#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                         NF-bcl2fastq
========================================================================================
 NF-bcl2fastq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/gnetsanet/NF-bcl2fastq
 #### Authors
 Netsanet Gebremedhin gnetsanet <gnetsanet@gmail.com> - https://github.com/gnetsanet>
----------------------------------------------------------------------------------------
*/

import groovy.util.*
import java.io.File


// To-do:
// 1. Add help message

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
params.dryrun = false
params.testmode = false // Perform end-to-end test with a test dataset - take a fraction of the data ( two tiles) from any instrument run
params.skipExtractBarcodes = false // Do not run the ExtractIlluminaBarcodes piece. Should add more and clever ways of skipping a given step. Useful in development
params.skipCreateBarcodeSummary = false
params.skipKronaReports = false
params.skipCombinedKronaReports = false
params.skipIlluminate = false
params.skipRunKrakenClassifyBam  = false
params.skipMergingLibraryBams = false
params.skipLinkLibraryBams = false
params.run_uri = '' // Initialize run URI
params.admin = params['emails'].admin

//=============================================== Staging and Preprocessing ====================================================

metrics_name="${params.max_no_calls}nc_${params.min_mismatch_delta}mmd_${params.max_mismatches}mis_bc_metrics.txt"
runDir = MiscUtils.getRunDirectory(samplesheet)
barcodesOutputDir = MiscUtils.createOutputDir(samplesheet, workflow.resume)
(tile_xml, instrument) = MiscUtils.getTileXmlFile(samplesheet)
(tilesList, preProcessDir) = MiscUtils.retrieveListOfTiles(tile_xml, instrument, runDir)
(headerLines, experimentName) = MiscUtils.parseSampleSheet(samplesheet)
(num_reads,run_id, flowcell, machine_name, lanecount, read_structure) = MiscUtils.parseRunXml(runDir, samplesheet, headerLines)


// Creates lane*multiplex_params.txt and lane*library_params.txt for each lane for use by picard. Will be replaced by fgbio's functions
MiscUtils.retrieveLibraryParams(preProcessDir, samplesheet, lanecount, headerLines, num_reads)
readGroups = MiscUtils.getReadGroups(samplesheet, headerLines, num_reads, lanecount)

// Convert text files to channels so they can be used as inputs in downstream processes
barcodeFiles = Channel.fromPath("${preProcessDir}/lane*_barcode_params.txt")

Channel.from(1..lanecount.toInteger())
	.into { lanes1; lanes2; lanes3; lanes4; lanes5; lanes6; lanes7; lanes8 }


Channel.fromPath(samplesheet)
  .splitCsv(header: true, skip: headerLines, strip: true)
	.map { row ->
		def m = row =~ /[a-zA-Z0-9_.]+\@[a-zA-Z]+.com/;
		m[0]
	}
	.unique()
	.into { user_emails1; user_emails2; user_emails3 }

// Preparing data structure for adding read group information to BAMs followed by merging into library BAMs.


readGroups
        .spread(lanes1)
        .map { tuple(it[0], it[1], it[2]) }
        .into { libraryParamsInfo1; libraryParamsInfo2; libraryParamsInfo3 }
        

def paired_end = (MiscUtils.getReadsToTransfer(num_reads).size()>1) ? true : false


// For merging library BAMS into a single lane BAM

libraries=Channel.fromPath("${preProcessDir}/lane${lanecount}_library_params.txt")
        .splitCsv(header: true, sep: '\t') 
        .map { line ->
        	line.SAMPLE_ALIAS
        }

lanes5 //Spread lanes over libraries
    .combine(libraries)
    .groupTuple()
    .set { lanesAndSamples }

// List of Libraries for Kraken Classification
Channel.fromPath("${preProcessDir}/lane${lanecount}_library_params.txt")
        .splitCsv(header: true, sep: '\t') 
        .map { line ->
        		line.OUTPUT
        }
        .into { sampleNamesForKrona; sampleNamesForKraken; sampleNamesForLink}


// For linking per lane combined BAMs and library BAMS
lanes7
    .spread(user_emails1)
    .set { laneAndUser }

sampleNamesForLink
	.spread(lanes8)
	.spread(user_emails3)
	.set { libraryLaneUser }

//End of Staging and Preprocessing =======================================================================================================

// Start of Processes ####################################################################################################################


process ExtractIlluminaBarcodes {

	tag { ["ExtractIlluminaBarcodes", lane ] }

	publishDir "${barcodesOutputDir}", mode: 'move', overwrite: false

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun

	input:
	file barcodeFile from barcodeFiles
	val lane from lanes4

	output:
	file ("*_barcode.txt") into ( barcodeTextFilesForSummary, barcodeTextFilesForBcl2Fastq, barcodeTextFilesForUnalignedBam)

	script:
	"""
	${params.java}/java -Xmx64000000k -jar ${params.picard}/picard.jar ExtractIlluminaBarcodes \
	MAX_NO_CALLS=${params.max_no_calls} MIN_MISMATCH_DELTA=${params.min_mismatch_delta} \
	MAX_MISMATCHES=${params.max_mismatches} NUM_PROCESSORS=${params.extract_bc_cpus} \
	read_structure=${read_structure} \
	LANE=${lane} \
	BASECALLS_DIR="${runDir}/Data/Intensities/BaseCalls/" \
	METRICS_FILE="${barcodesOutputDir}/L${lane}_${metrics_name}" \
	BARCODE_FILE=${barcodeFile} \
	OUTPUT_DIR=.
	"""
}


process createBarcodeSummary {

	errorStrategy 'retry'
	maxRetries 3

	when: !params.dryrun

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


process produceUnalignedBam {

	tag { [lane] }

	publishDir "${barcodesOutputDir}/UnalignedBAM/L${lane}/", mode: 'move', overwrite: false

	errorStrategy 'retry'
	maxRetries 3

	echo true

	when: !params.dryrun

	input:
	val lane from lanes6
	file barcodeText from barcodeTextFilesForUnalignedBam.toList()

	output:
	file("*.bam") into rawBams


	script:
	"""
	${params.java}/java -Xmx16000000k -jar ${params.picard}/picard.jar IlluminaBasecallsToSam \
	NUM_PROCESSORS=${params.tile_process_cpus} \
	read_structure=$read_structure \
	RUN_BARCODE=${machine_name} \
	LANE=${lane} \
	BASECALLS_DIR="${runDir}/Data/Intensities/BaseCalls/" \
	BARCODES_DIR=${barcodesOutputDir} \
	LIBRARY_PARAMS="${preProcessDir}/lane${lanecount}_library_params.txt"
	"""
}


process addReadGroups {

	tag { [library,  lane] }

	publishDir "${barcodesOutputDir}/CombinedBAMs", mode: 'move', overwrite: false

	errorStrategy 'retry'
	maxRetries 3

	echo true

	when: !params.dryrun

	input:
	set library, rdg, lane from libraryParamsInfo1
	file("*.bam") from rawBams.toList()

	output:
	file("out_L${lane}_${library}.bam") into (unalignedBamWithRG )


	script:
	"""
	${params.java}/java -Xmx16000000k -jar ${params.picard}/picard.jar AddOrReplaceReadGroups \
	INPUT=${barcodesOutputDir}/UnalignedBAM/L${lane}/L${lane}_${library}.bam \
	OUTPUT=out_L${lane}_${library}.bam \
	RGID=${rdg} \
	RGLB=${rdg} \
	RGPL="illumina" \
	RGPU=${rdg} \
	RGSM=${library}
	"""
}


process renameBams {
	tag { [library,  lane] }

	publishDir "${barcodesOutputDir}/CombinedBAMs", mode: 'move', overwrite: false

	errorStrategy 'retry'
	maxRetries 3

	echo true

	when: !params.dryrun

	input:
	set library, rdg, lane from libraryParamsInfo2
	file(rdgLibraryBam) from unalignedBamWithRG.collect()

	output:
	file("L${lane}_${library}.bam") into ( renamedBams1, renamedBams2, renamedBams3)

	script:
	"""
	cp ${barcodesOutputDir}/CombinedBAMs/out_L${lane}_${library}.bam L${lane}_${library}.bam
	rm ${barcodesOutputDir}/CombinedBAMs/out_L${lane}_${library}.bam
	"""
}


process mergeSampleBamsIntoLane {
	tag { [lane] }

	publishDir "${barcodesOutputDir}/CombinedBAMs/", mode: 'move', overwrite: false

	errorStrategy 'retry'
	maxRetries 5
	echo true

	when: !params.dryrun && !params.skipMergingLibraryBams

	input:
	set lane, libraries from lanesAndSamples
	file(libraryBam) from renamedBams1.collect()

	output:
	file("L${lane}_combined.bam") into laneCombinedBam

	script:
	inputBams = "INPUT=${barcodesOutputDir}/CombinedBAMs/" +  libraryBam.join(" INPUT=${barcodesOutputDir}/CombinedBAMs/")
	"""
	${params.java}/java -Xmx16000000k -jar ${params.picard}/picard.jar MergeSamFiles \
	${inputBams} \
	OUTPUT=L${lane}_combined.bam
	"""
}


process linkLibraryBams {

	tag {[library, lane, user ]}

	publishDir "${barcodesOutputDir}", mode: 'move', overwrite: false

	errorStrategy 'retry'
	maxRetries 5

	when: !params.dryrun && !params.skipLinkLibraryBams // we do not need to link per library bam for the purity assessment only combined BAM, so make this process optional

	input:
	set library, lane, user from libraryLaneUser
	file unAlnBam from renamedBams2.toList()

	output:
	file("${library}_${user}_Link_Complete.txt") into bamLinkCompletionMarkers1

	script:
	"""
	mkdir -p ${params.symlink_dest}/${user}/${run_id}_shep && chmod 777 ${params.symlink_dest}/${user}/${run_id}_shep
	ln -sf  ${barcodesOutputDir}/CombinedBAMs/${library} ${params.symlink_dest}/${user}/${run_id}_shep
	touch ${library}_${user}_Link_Complete.txt
	"""

}


process linkLaneBams {

	tag {[lane, user]}

	publishDir "${barcodesOutputDir}/Flags", mode: 'move', overwrite: false

	errorStrategy 'retry'
	maxRetries 5

	when: !params.dryrun && !params.skipMergingLibraryBams

	input:
	set lane, user from laneAndUser
	file laneBam from laneCombinedBam

	output:
	file("${lane}_Link_Complete.txt") into bamLinkCompletionMarkers2

	script:
	if('miseq'.equalsIgnoreCase(instrument)) { fcell = flowcell.split('-')[1] } else { fcell = flowcell } // Name the combined BAM after flowcell
	
	"""
	mkdir -p ${params.symlink_dest}/${user}/${run_id}_shep && chmod 777 ${params.symlink_dest}/${user}/${run_id}_shep
	ln -sf  ${barcodesOutputDir}/CombinedBAMs/L${lane}_combined.bam ${params.symlink_dest}/${user}/${run_id}_shep/${fcell}.bam
	touch ${lane}_Link_Complete.txt
	"""

}

process runKrakenClassifyBam {
	// Run kraken classify using unaligned BAMs as inputs
	tag { sampleId }

	publishDir "${barcodesOutputDir}/Flags" , mode: 'move', overwrite: false

	errorStrategy 'retry'
	maxRetries 5

	when: !params.dryrun && !params.skipRunKrakenClassifyBam

	input:
	val sampleId from sampleNamesForKraken
	file unAlnBam from renamedBams3.toList()

	output:
	file("${sampleId}_Kraken_Classify_Complete.txt") into ( krakenClassifyBamCompletionMarkerKrona, krakenClassifyBamCompletionMarkerKronaCombined )

	script:
		if(paired_end) {
						readNum = MiscUtils.getReadsToTransfer(num_reads)[-1]
						"""
						mkdir -p ${barcodesOutputDir}/Kraken
						samtools sort -n ${barcodesOutputDir}/CombinedBAMs/${sampleId} -o ${barcodesOutputDir}/CombinedBAMs/${sampleId}.sorted
						bedtools bamtofastq -i ${barcodesOutputDir}/CombinedBAMs/${sampleId}.sorted -fq ${barcodesOutputDir}/CombinedBAMs/${sampleId}.1.fastq -fq2 ${barcodesOutputDir}/CombinedBAMs/${sampleId}.${readNum}.fastq
						ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no seq-shepherd@seq-himem02 '${params.kraken} --db  ${params.krakenDB} --paired ${barcodesOutputDir}/CombinedBAMs/${sampleId}.1.fastq ${barcodesOutputDir}/CombinedBAMs/${sampleId}.${readNum}.fastq --threads 4 --output ${barcodesOutputDir}/Kraken/${sampleId}.kraken'
						touch ${sampleId}_Kraken_Classify_Complete.txt
						"""
		} else {
						"""
						mkdir -p ${barcodesOutputDir}/Kraken
						samtools sort -n ${barcodesOutputDir}/CombinedBAMs/${sampleId} -o ${barcodesOutputDir}/CombinedBAMs/${sampleId}.sorted
						bedtools bamtofastq -i ${barcodesOutputDir}/CombinedBAMs/${sampleId}.sorted -fq ${barcodesOutputDir}/CombinedBAMs/${sampleId}.1.fastq
						ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no seq-shepherd@seq-himem02 '${params.kraken} --db ${params.krakenDB} ${barcodesOutputDir}/CombinedBAMs/${sampleId}.1.fastq --threads 4 --output ${barcodesOutputDir}/Kraken/${sampleId}.kraken'
						touch ${sampleId}_Kraken_Classify_Complete.txt
						"""
		}

}


process createkronaReports {

	tag { "${sampleId}" }

	publishDir "${barcodesOutputDir}/Flags" , mode: 'move', overwrite: false


	errorStrategy 'retry'
  	maxRetries 5

	when: !params.dryrun && !params.skipRunKrakenClassifyBam


    input:
	val sampleId from sampleNamesForKrona
	file krakenMarker from krakenClassifyBamCompletionMarkerKrona.toList()


	output:
	file("*_kronaReport_Complete.txt") into kronaReport

    script:
	"""
	mkdir -p ${barcodesOutputDir}/Krona
	ssh  -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no seq-shepherd@seq-himem02 '${params.ktImportTaxonomy} -o ${barcodesOutputDir}/Krona/${sampleId}.krona.html  -t 3 -s 4 ${barcodesOutputDir}/Kraken/${sampleId}.kraken'
	touch ${sampleId}_kronaReport_Complete.txt
    """

}


process createCombinedKronaReport {

	tag { "Creating a combined/per run Kraken Report."}

 	publishDir "${barcodesOutputDir}/Flags" , mode: 'move', overwrite: false

	errorStrategy 'retry'
	maxRetries 5

	when: !params.dryrun && !params.skipRunKrakenClassifyBam

	input:
	file krakenMarkerCombined from krakenClassifyBamCompletionMarkerKronaCombined.toList()

	output:
	file("Combined_Krona_Report_Complete.txt") into kronaCombinedOut

	script:
	krona_inputs = Channel.fromPath("${barcodesOutputDir}/Kraken/**.kraken").toList().getVal().join(" ")
	if (!params.skipRunKrakenClassifyBam) {
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

process extractRunMetrics {

	tag { 'Parsing InterOp/Illuminate'}

	publishDir "${barcodesOutputDir}", mode: 'move', overwrite: false

	errorStrategy 'retry'
	maxRetries 5

	when: !params.dryrun

	output:
	file('runMetrics.txt') into runMetrics

	script:
	"""
	illuminate --tile --quality ${runDir} > runMetrics.txt
	"""
}


process notifyUsersRuby {

		tag { userEmailR }
		publishDir "${barcodesOutputDir}/Flags", mode: 'move', overwrite: false

		errorStrategy 'retry'
   		maxRetries 5
   		echo true

		when: !params.nomail && !params.dryrun

		input:
		val userEmailR from user_emails2
		file kronaR from kronaCombinedOut
		file illuminate from runMetrics

		script:
		def filesToAttach = []
		if(!params.skipMultiQc) { filesToAttach.add("${barcodesOutputDir}/MultiQC/multiqc_report.html")	}									
		if(!params.skipRunKrakenClassifyBam ) { filesToAttach.add("${barcodesOutputDir}/Contamination.kraken_Report.html")}  
		attachments = filesToAttach.join(", ")
		// Barcode frequency will always be calculated
		// However, MulitQC and Combined Kraken/Krona Report
		// generation may not/should not always happen
		// Hence, those MultiQC and Combined Krona files should be added only if the processes generating them have not been skipped
		"""
		ruby ${baseDir}/../create_email_message.rb -b ${barcodesOutputDir}/L*_${metrics_name} \
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