#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
						 NF-bcl2fastq
========================================================================================
 NF-bcl2fastq Analysis Pipeline.
----------------------------------------------------------------------------------------
*/

import groovy.util.*
import java.io.File


// To-do:
// 1. Java memory allocation is hard-coded right now à la -Xmx64000000k - Not super urgent since we can specify resources per process.
// 2. Add help message


// Defaults and input params ****************************************************************************************************
samplesheet = params.samplesheet
params.max_mismatches = params['cutoffs_and_thresholds'].max_mismatches
params.max_no_calls = params['cutoffs_and_thresholds'].max_no_calls
params.min_mismatch_delta = params['cutoffs_and_thresholds'].min_mismatch_delta
params.extract_bc_cpus = params['cutoffs_and_thresholds'].extract_bc_cpus
params.tile_process_cpus = params['cutoffs_and_thresholds'].tile_process_cpus
params.picard = params['script_paths'].picard
params.java = params['script_paths'].java
params.post_run_scripts = params['script_paths'].post_run_scripts
params.fasqtc = params['script_paths'].fastqc
params.multiqc = params['script_paths'].multiqc
params.kraken = params['script_paths'].kraken
params.krakenDB = params['script_paths'].krakenDB
params.ktImportTaxonomy = params['script_paths'].ktImportTaxonomy
params.symlink_dest = params['other_paths'].symlink_dest
params.link = true // whether or not to create symlinks to combined fastqs. Set to false while testing since seq-shepherd, too does this.
params.nomail = false // Good parameter to have for testing since we do not want to spam people
params.mail_function = 'ruby'
params.dryrun = false
params.testmode = false // Perform end-to-end test with a test dataset - take a fraction of the data ( two tiles) from any instrument run
params.skipExtractBarcodes = false // Do not run the ExtractIlluminaBarcodes piece. Should add more and clever ways of skipping a given step. Useful in development
params.skipBcl2Fastq = false
params.skipCopyCombineFastq = false
params.skipCreateBarcodeSummary = false
params.skipKrakenClassify = false
params.skipKronaReports = false
params.skipCombinedKronaReports = false
params.skipMultiQc = false
params.skipIlluminate = false
params.run_uri = '' // Initialize run URI
params.admin = params['emails'].admin

//=============================================== Staging and Preprocessing ====================================================

metrics_name="${params.max_no_calls}nc_${params.min_mismatch_delta}mmd_${params.max_mismatches}mis_bc_metrics.txt" //Construct name for the metrics file
runDir = MiscUtils.getRunDirectory(samplesheet) // Parse run directory base name from samplesheet path
barcodesOutputDir = MiscUtils.createOutputDir(samplesheet, workflow.resume) // Output directory where final and intermediate files are stored
(tile_xml, instrument) = MiscUtils.getTileXmlFile(samplesheet) // Return xml file, location and name of xml could be different in different instruments
(tilesList, preProcessDir) = MiscUtils.retrieveListOfTiles(tile_xml, instrument, runDir) //List of tiles to parallelize processes on/Nextflow Channel 
(headerLines, experimentName) = MiscUtils.parseSampleSheet(samplesheet)
(num_reads,run_id, flowcell, machine_name, lanecount, read_structure, has_umi) = MiscUtils.parseRunXml(runDir, samplesheet, headerLines) // Gather basic run data.



// Creates lane*multiplex_params.txt and lane*library_params.txt for each lane for use by picard. Will be replaced by fgbio's functions
MiscUtils.retrieveLibraryParams(preProcessDir, samplesheet, lanecount, headerLines, num_reads)

// Convert text files to channels so they can be used as inputs in downstream processes
barcodeFiles = Channel.fromPath("${preProcessDir}/lane*_barcode_params.txt")

// Create an input that can be consumed by 'copyCombineFastq' - a combination of sample and read { 1, 2} to combined per tile fastqs into

Channel.fromPath(samplesheet)
	.splitCsv(header: true, skip: headerLines, strip: true)
	.map { row ->
		libraryName = MiscUtils.reformatSampleId(row.Sample_ID) //Replace one or more spaces and backticks within library names
	}
	.unique()
	.into { sampleNames; sampleNamesForKraken; sampleNamesForKrona }

sampleNames
	.spread(Channel.from(MiscUtils.getReadsToTransfer(num_reads, has_umi)))
	.into { reads2Transfer; readsForFastqc}

Channel.from(1..lanecount.toInteger())
	.into { lanes; lanesForMapping; lanesForMergingUnBams}

// Create a list of  available lane and tile combinations that can be used to parallelize the bcl2fastq process
// tilesList is returned by a function that parses the run xml file
Channel.from(tilesList)
	.splitCsv(header: false, limit: params.testmode ? 3 : 0) // Counter-intuitively, limit of 0 means everything/No limit.
	.map {
		lane= (instrument == 'novaseq') ? it[0].split("_")[0] : 1
		tile= (instrument == 'novaseq') ? it[0].split("_")[1] : it[0]
		[lane, tile]
	}
	.set { tilesFastq }

// Two channel subsetting and transformation to generate inputs for use by linkFastqsToUsersFolders
// and notifyUsers processes

Channel.fromPath(samplesheet)
		.splitCsv(header: true, skip: headerLines, strip: true)
		.map { row ->
				sampleId = MiscUtils.reformatSampleId(row.Sample_ID) //Strip a csv field's leading and trailing white spaces with splitCsv's 'strip : true' option. Also replace backticks with -, seen in 180426_M05473_0127_000000000-BRYR4/SampleSheet.csv
				def m = row =~ /[a-zA-Z0-9_.]+\@[a-zA-Z]+.com/; // project owner's email may not be in a specific column albeit it usually is in a column called 'Description'
								user = m[0]
				[sampleId, user]
		}
		.into { sampleAndUser; sampleAndUser2; sampleAndUser3; sampleAndUser4 }

Channel.fromPath(samplesheet)
	.splitCsv(header: true, skip: headerLines, strip: true)
	.map { row ->
		def m = row =~ /[a-zA-Z0-9_.]+\@[a-zA-Z]+.com/;
		m[0]
	}
	.unique()
	.into { user_emails1; user_emails2 }

Channel
	.from(file("${preProcessDir}/lane${lanecount}_barcode_params.txt").readLines().drop(1))
	.map { line ->
					cols = line.split("\t")
					cols[-1] // There may not always be a barcode name in the lane1_barcode_params.txt so better to use the last column to get to library names
	}
	.spread(lanesForMergingUnBams)
	.spread(file("${preProcessDir}/tilesList.csv").readLines())
	// .subscribe { println "$it" }
	.into { mappingChannel; mappingChannelForMergingBams }

mappingChannelForMergingBams
		.map { tuple(it[0], it[1] + "_" + it[2]) }
		.groupTuple()
		//.subscribe { println it }
		.set { tileIdsMergeBams }

// Preparing input for merging unaligned BAMs
Channel
	.from(file("${preProcessDir}/lane${lanecount}_library_params.txt").readLines().drop(1))
	.map { line ->
					cols = line.split("\t")
					cols[-1]
	}
	.spread(lanesForMapping)
	.spread(file("${preProcessDir}/tilesList.csv").readLines())
	.map { tuple(it[0], it[1] + "_" + it[2]) }
	.groupTuple()
	.set { tileIdsMergeUnalignedBams }

def paired_end = (MiscUtils.getReadsToTransfer(num_reads, has_umi).size()>1) ? true : false

//End of Staging and Preprocessing =======================================================================================================

// Start of Processes ####################################################################################################################


process ExtractIlluminaBarcodes {

	tag { ["ExtractIlluminaBarcodes", lane ] }

	publishDir "${barcodesOutputDir}", mode: 'copy', overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun && !params.skipExtractBarcodes

	input:
	file barcodeFile from barcodeFiles
	val lane from lanes

	output:
	file ("extracting_lane${lane}_complete.txt") into ( barcodeTextFilesForSummary, barcodeTextFilesForBcl2Fastq)

	script:
	"""
	${params.java}/java -Xmx64000000k -jar ${params.picard}/picard.jar ExtractIlluminaBarcodes \
	MAX_NO_CALLS=${params.max_no_calls} MIN_MISMATCH_DELTA=${params.min_mismatch_delta} \
	MAX_MISMATCHES=${params.max_mismatches} NUM_PROCESSORS=${params.extract_bc_cpus} \
	read_structure=${read_structure} \
	LANE=${lane} \
	BASECALLS_DIR="${runDir}/Data/Intensities/BaseCalls/" \
	METRICS_FILE="${barcodesOutputDir}/L${lane}_${metrics_name}" \
	BARCODE_FILE=${preProcessDir}/lane${lane}_barcode_params.txt
	touch extracting_lane${lane}_complete.txt
	"""
}


process createBarcodeSummary {

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun && !params.skipCreateBarcodeSummary

	input:
	file barcodeText from barcodeTextFilesForSummary.toList()

	script:
	Channel.fromPath("${barcodesOutputDir}/L*_${metrics_name}")
		.first()
		.collectFile(skip: 6) // ignore header which is the first 6 lines
		.splitText()
		.filter({it.tokenize('\t').size() >= 3}) // split tab-delimited file by tab and ignore blank lines using the count of column data as a proxy for blankness.
		.map { row ->
			columns=row.tokenize("\t")

			if (columns[0] =~ /^N+$/) {
				columns[2] = 'Unknown' // rename columns since no library or sample name for unassigned barcodes and none is output in the file we are processing
				columns[3] = 'Undetermined'
			}
			[columns[0], columns[1], columns[2], columns[3], columns[4], columns[11]].join("\t")

		}
		.collectFile(
			name: "${barcodesOutputDir}/barcode_summary.txt", newLine: true, sort: false
		)
		'''
		'''
}


process bcl2fastq {

	tag { [lane, tile] }

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun

	input:
	set lane, tile from tilesFastq
	file barcodeText from barcodeTextFilesForBcl2Fastq.collect()

	output:
	file("*.fastq.gz") into fastqFiles


	script:
	"""
	${params.java}/java -Xmx16000000k -jar ${params.picard}/picard.jar IlluminaBasecallsToFastq \
	NUM_PROCESSORS=${params.tile_process_cpus} \
	read_structure=$read_structure \
	RUN_BARCODE=${machine_name} \
	LANE=${lane} \
	FIRST_TILE=${tile} \
	TILE_LIMIT=1 \
	MACHINE_NAME=${machine_name} \
	FLOWCELL_BARCODE=${flowcell} \
	BASECALLS_DIR="${runDir}/Data/Intensities/BaseCalls/" \
	MULTIPLEX_PARAMS="${preProcessDir}/lane${lane}_multiplex_params.txt" \
	MAX_READS_IN_RAM_PER_TILE=3000000 \
	MAX_RECORDS_IN_RAM=3000000 \
	COMPRESS_OUTPUTS=true \
	COMPRESSION_LEVEL=1 \
	TMP_DIR=/mnt/flash_scratch/tmp
	"""
}


fastqFiles
	.toList().flatten()
	.map {
		file -> tuple(file.baseName.replaceAll(/^L\d+_/,''), file)
	}
	.groupTuple()
	.set { fastqFilesForCombining }


process copyCombineFastq {

	tag { libraryRead }
	publishDir "${barcodesOutputDir}/combined_fastq" , mode: 'copy', overwrite: true

	errorStrategy 'retry'
	maxRetries 1
	echo true

	when: !params.dryrun

	input:
	set val(libraryRead), val(fastqs) from fastqFilesForCombining

	output:
	file("*.fastq.gz") into ( combineFastqCompletionMarkers, combineFastqCompletionMarkersForKraken)

	script:
	def libraryName = libraryRead.replaceAll(/^L\d+_/,'')
	"""
	zcat ${fastqs.join(" ")}| /mnt/galaxy/data/galaxy/sw/bin/pigz -p 8 > ${libraryName}.gz
	"""
}


process runKrakenClassify {

	tag { sampleId }

	publishDir "${barcodesOutputDir}/Flags" , mode: 'copy', overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun

	input:
	val sampleId from sampleNamesForKraken
	file combinecompletionMarker from combineFastqCompletionMarkersForKraken.toList()

	output:
	file("*_Kraken_Classify_Complete.txt") into ( krakenClassifyCompletionMarkerKrona, krakenClassifyCompletionMarkerKronaCombined )

	script:
	if(!params.skipKrakenClassify) {
			if(paired_end) {
							readNum = MiscUtils.getReadsToTransfer(num_reads, has_umi)[-1]
					"""
					mkdir -p ${barcodesOutputDir}/Kraken
					ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no seq-shepherd@seq-himem02 << EOF
					exec 200>/tmp/NF_krakenClassify.lock || exit 1
					flock 200 || exit 1
					${params.kraken} --db  ${params.krakenDB} --paired ${barcodesOutputDir}/combined_fastq/${sampleId}.1.fastq.gz ${barcodesOutputDir}/combined_fastq/${sampleId}.${readNum}.fastq.gz --threads 4 --output ${barcodesOutputDir}/Kraken/${sampleId}.kraken
					flock -u 200
					EOF
					touch ${sampleId}_Kraken_Classify_Complete.txt
					"""
			} else {
					// It is better to create a new matrix spread that contains combination of library id and fastq paths than repeat ${barcodesOutputDir}/combined_fastq/
					"""
					mkdir -p ${barcodesOutputDir}/Kraken
					ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no seq-shepherd@seq-himem02 << EOF
					exec 200>/tmp/NF_krakenClassify.lock || exit 1
					flock 200 || exit 1
					${params.kraken} --db ${params.krakenDB} ${barcodesOutputDir}/combined_fastq/${sampleId}.*.fastq.gz --threads 4 --output ${barcodesOutputDir}/Kraken/${sampleId}.kraken
					flock -u 200
					EOF
					touch ${sampleId}_Kraken_Classify_Complete.txt
					"""
			}
		} else {
			"""
			touch ${sampleId}_Kraken_Classify_Complete.txt
			"""
		} 	
}


process createkronaReports {

	tag { "Creating Krona html reports for library ${sampleId}" }

	publishDir "${barcodesOutputDir}/Flags" , mode: 'copy', overwrite: true


	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun && !params.skipKrakenClassify && !params.skipKronaReports


	input:
	val sampleId from sampleNamesForKrona
	file krakenMarker from krakenClassifyCompletionMarkerKrona.toList()


	output:
	file("*_kronaReport_Complete.txt") into kronaReport

	script:
	"""
	mkdir -p ${barcodesOutputDir}/Krona
	ssh  -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no seq-shepherd@seq-himem02 << EOF
	exec 200>/tmp/NF_kronaReport.lock || exit 1
	flock 200 || exit 1
	${params.ktImportTaxonomy} -o ${barcodesOutputDir}/Krona/${sampleId}.krona.html  -t 3 -s 4 ${barcodesOutputDir}/Kraken/${sampleId}.kraken
	flock -u 200
	EOF
	touch ${sampleId}_kronaReport_Complete.txt
	"""

}


process zipKronaReports {

	tag { "Archiving Krona Reports"}

	publishDir "${barcodesOutputDir}/Krona" , mode: 'copy', overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun && !params.skipKrakenClassify && !params.skipKronaReports


	input:
	file(kronaEndMarker) from kronaReport.toList()

	output:
	file("${run_id}.krona.zip") into zipOut

	script:
	krona_reports = Channel.fromPath("${barcodesOutputDir}/Krona/**.html").toList().getVal().join(" ")
	"""
	zip -j  ${run_id}.krona.zip  ${krona_reports}
	"""

}


process createCombinedKronaReport {

	tag { "Creating a combined/per run Kraken Report."}

	publishDir "${barcodesOutputDir}/Flags" , mode: 'copy', overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun 

	input:
	file krakenMarkerCombined from krakenClassifyCompletionMarkerKronaCombined.toList()

	output:
	file("Combined_Krona_Report_Complete.txt") into ( kronaCombinedOut1, kronaCombinedOut2 )

	script:
	krona_inputs = Channel.fromPath("${barcodesOutputDir}/Kraken/**.kraken").toList().getVal().join(" ")
	if (!params.skipCombinedKronaReports && !params.skipKrakenClassify) {
			"""
			ssh seq-shepherd@seq-himem02 '${params.ktImportTaxonomy} -o ${barcodesOutputDir}/Contamination.kraken_Report.html -t 3 -s 4  ${krona_inputs}'
			touch Combined_Krona_Report_Complete.txt
			"""
	} else {
		"""
		touch 'Combined_Krona_Report_Complete.txt' # Create a bogus completion flag file so that subsequent processes that depend on this process ( if any) get executed.
		"""
	}

}


process runFastQc {

	tag { sampleId }

	publishDir "${barcodesOutputDir}/Flags" , mode: 'copy', overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun

	input:
	set sampleId, user from sampleAndUser2
	file combinedFastqs from combineFastqCompletionMarkers


	output:
	file("${sampleId}_FastQC_Complete.txt") into fastQCreports

	script:
	"""
	mkdir -p ${barcodesOutputDir}/fastqc
	${params.fasqtc} --extract  -o ${barcodesOutputDir}/fastqc -f fastq -q ${barcodesOutputDir}/combined_fastq/${sampleId}.*.fastq.gz
	touch ${sampleId}_FastQC_Complete.txt
	"""
}


process runMultiQc {

	tag { "multiQC"}

	publishDir "${barcodesOutputDir}/MultiQC", mode: 'copy', overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun && !params.skipMultiQc

	input:
	file fastqcReport from fastQCreports.toList()

	output:
	file("multiqc_report.html") into multiQcReport

	script:
	"""
	${params.multiqc} ${barcodesOutputDir}/fastqc  -n  multiqc_report.html -o ./  --force
	"""
}


process linkFastqsToUsersFolders {

	tag {[user, sampleId]}

	publishDir "${barcodesOutputDir}/Flags", mode: 'copy', overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	when: params.link && !params.dryrun

	input:
	set sampleId, user from sampleAndUser
	file multiQc from multiQcReport
	file krona from kronaCombinedOut2

	output:
	file("${sampleId}_${user}_Link_Complete.txt") into linkCompletionMarkers

	script:

	"""
	mkdir -p ${params.symlink_dest}/${user}/${run_id}_shep && chmod 777 ${params.symlink_dest}/${user}/${run_id}_shep
	ln -sf ${barcodesOutputDir}/combined_fastq/${sampleId}*.fastq.gz ${params.symlink_dest}/${user}/${run_id}_shep
	touch ${sampleId}_${user}_Link_Complete.txt
	"""

}


process calculateBarcodeFrequency {

	tag { 'Calculating unassigned barcode frequency'}

	publishDir "${barcodesOutputDir}", mode: 'copy', overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun

	input:
	file linkCompletionMarker from linkCompletionMarkers.toList()

	output:
	file('unassigned_barcode_frequency.txt') into ( barcode_distribution1, barcode_distribution2, barcode_distribution3) // We need 2 copies of barcode_distribution

	script:
	"""
	# Do not make the assumption there will always be tile 1101 and lane 1, particularly not the lane.
	# For instance, in NextSeq tile patterns are more like 11101 as opposed to 1101 for novaseq since there are a different count of lanes, surface, swath, tile in NextSeq
	zcat ${barcodesOutputDir}/fastq/*/L*_unassigned.barcode_*.fastq.gz |head -n 100000|sed -n '2~4p'|sort|uniq -c|sort -rn|head>unassigned_barcode_frequency.txt
	"""
}


process extractRunMetrics {

	tag { 'Parsing InterOp/Illuminate'}

	publishDir "${barcodesOutputDir}", mode: 'copy', overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun

	output:
	file('runMetrics.txt') into runMetrics

	script:
	"""
	illuminate --tile --quality ${runDir} > runMetrics.txt
	"""
}


process prepareEmailMessage {

	tag { 'Preparing email messages'}

	publishDir "${barcodesOutputDir}", mode: 'copy', overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun

	input:
	file barcodeDistribution from barcode_distribution3.toList()

	output:
	file('message.txt') into (email_message1, email_message2)

	script:
	"""
	printf "Run ${experimentName} ( ${flowcell}) processing is complete. Output files have been linked to your user folder.\n\n \
	Output files have been linked to your user folder ${params.symlink_dest}/${user}/${run_id}_shep.\n\n \
	Barcode Summary from Lane 1 and the 10 most frequent barcodes observed in unassigned fastq file (first 100k reads) are attached.\n\n \
	Please remember to clean the instrument.">message.txt
	"""

}


process notifyUsersRuby {

	tag { userEmailR }
	publishDir "${barcodesOutputDir}/Flags", mode: 'copy', overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	when: !params.nomail && params.mail_function.toLowerCase()=="ruby" && !params.dryrun

	input:
	val userEmailR from user_emails2
	file barcodeFreqR from barcode_distribution2
	file messageR from email_message2
	file kronaR from kronaCombinedOut1
	file illuminate from runMetrics
	file zip from zipOut

	script:
	// def filesToAttach = ["${barcodesOutputDir}/unassigned_barcode_frequency.txt"] 
	def filesToAttach = []
	if(!params.skipMultiQc) { filesToAttach.add("${barcodesOutputDir}/MultiQC/multiqc_report.html")	}									
	if(!params.skipCombinedKronaReports && !params.skipKrakenClassify ) {
	filesToAttach.add("${barcodesOutputDir}/Contamination.kraken_Report.html")
	}  
	if(!params.zipKronaReports && !params.skipKrakenClassify) {
		filesToAttach.add("${barcodesOutputDir}/Krona/${run_id}.krona.zip")
	}
	attachments = filesToAttach.join(", ")
	// Barcode frequency will always be calculated
	// However, MulitQC and Combined Kraken/Krona Report
	// generation may not/should not always happen
	// Hence, those MultiQC and Combined Krona files should be added only if the processes generating them have not been skipped
	"""
	ruby ${baseDir}/../create_email_message.rb -b ${barcodesOutputDir}/L*_${metrics_name} \
	-u ${barcodesOutputDir}/unassigned_barcode_frequency.txt \
	-i ${barcodesOutputDir}/runMetrics.txt \
	-e \"${experimentName}\" \
	-a ${userEmailR} \
	-f ${run_id} \
	-l ${params.run_uri}
	ruby ${baseDir}/../mail.rb -r "${userEmailR},${params.admin}" \
	-s  "Run Complete: ${experimentName} (${run_id})" -b message.html -f \"${attachments}\"
	touch 'pipeline_complete.txt' # Marks the completion of the pipeline processs
	"""
}


workflow.onError {
	
	def subject = "Seq-shepherd Failure: ${run_id}"
	def recipients = params.admin


	['mail', '-s', subject, recipients].execute() << """

	Pipeline execution summary
	---------------------------
	Completed at : ${workflow.complete}
	Duration     : ${workflow.duration}
	Success      : ${workflow.success}
	workDir      : ${workflow.workDir}
	exit status  : ${workflow.exitStatus}
	Error report : ${workflow.errorReport ?: '-'}
	Error Message: ${workflow.errorMessage ?: '-'}
	Commands     : ${workflow.commandLine}
	Run name     : ${workflow.runName}	 		
	"""

}
