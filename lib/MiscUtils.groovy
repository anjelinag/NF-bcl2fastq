import static nextflow.Nextflow.file
import nextflow.Channel


class MiscUtils {

	static def reformatSampleId(sampleId) {
		return sampleId.replaceAll("\\s","-").replaceAll(/[^A-Za-z0-9._-]/,"-") // Remove illegal characters since users use characters that are impermissible for file names eg. ampersand
	}

	static def getReadsToTransfer(num_reads, has_umi) {
	    def readsToTransfer
	    switch (num_reads) {
	        case 1:
	            readsToTransfer = [1]
	            break
	        case 2:
	            readsToTransfer = [1]
	            break
	        case 3:
	            readsToTransfer = [1,2] //The second read could be the third to be sequenced but still it is the second of two/paired end reads
	            break
	        case 4: // Dual indexed runs where there are a total of 4 reads are tricky as they may or may not have UMIs
	            if(has_umi) { 
	            	readsToTransfer = [1,2,3] // We need both read 1, UID and read 3 in runs containing dual barcodes and has UMIs
	        	} else {
	        		readsToTransfer = [1,2] // if the run has no UMI's, we do not need read 2.
	        	}
	            break
	    }
	    return readsToTransfer
	}

	static def getRunDirectory(filePath) {
		def x = filePath.tokenize("/")
		x.pop()
		return '/' +  x.join("/")

	}

	static def createOutputDir(filePath, resume) {

		def barcodesDir = new File(getRunDirectory(filePath) + "/Nextflow/Barcodes")
		if( resume ) { return barcodesDir  }
		if(!barcodesDir.exists()) {
			barcodesDir.mkdirs()
		} else {
			// Nextflow must have run before, so rename the outer directory
			// to-do: clean-up repetitive use of getRunDirectory
			new File(getRunDirectory(filePath) + "/Nextflow").renameTo(getRunDirectory(filePath) + "/Nextflow" + '_' + randomGenerator(7)) // rename the original outer 'Nextflow' directory
			barcodesDir.mkdirs()
		}
		return barcodesDir

	}

	static def getTileXmlFile(filePath) {
		// A function that is meant to identify the xml file that contains list of tiles. Tiles are listed in RunInfo.xml in the
		// case of Novaseq and NextSeq while in config.xml in the case of MiSeqs.
		def x = filePath.tokenize("/")
		x.pop()
		def run_directory = x.join("/")


		def runInfoXml = "/" + run_directory + "/RunInfo.xml"
		def configXml = "/" + run_directory + "/Data/Intensities/BaseCalls/config.xml" // For Miseqs the tile info is in config.xml file buried under the Basecalls dir.

		if ( file(configXml).exists() ) { return [configXml, 'miseq'] } // Although both Novaseq and miseq have RunInfo.xml, this would suffice to tell which
		if ( file(runInfoXml).exists() ) { return [runInfoXml, 'novaseq'] } // This could also be a nextseq, however miseq is the only one with a different xml file content

		exit 1, "Missing the expected run xml file: filePath"
	}

	static def parseRunXml(run_directory, ss, header) {
		def xmlFilePath = "/" + run_directory + "/RunInfo.xml" // The content varies but both Novaseq and Miseq have RunInfo.xml.
		def runInfo = new XmlParser().parse(new File(xmlFilePath))
		def num_reads = runInfo.Run.Reads.Read.size()
		def run_id = runInfo.Run.@'Id'[0]
		def flowcell = runInfo.Run.Flowcell.text()
		def machine_name = runInfo.Run.Instrument.text()
		def lanecount = runInfo.Run.FlowcellLayout.@'LaneCount'[0]

		def ssFile = new File(ss)
		def ssLines = ssFile.readLines()
		def ssColumnNames = ssLines.get(header.toInteger()) // Extract the column headers
		def index = getBarcodeColumns(ssColumnNames,'index') // Get the barcode/index column given all the column names in the  Sample Sheet
		def index2 = getBarcodeColumns(ssColumnNames,'index2')
		def has_umi = false

		// read_structure = ""
		// runInfo.Run.Reads.Read.each { r ->
		// 	if ('Y'.equalsIgnoreCase(r.@'IsIndexedRead')) {
		// 		read_structure += r.@'NumCycles' + "B"
		// 	} else {
		// 		read_structure += r.@'NumCycles' + "T"
		// 	}
		//
		// }
		def read_structure = []

		// Gather the templates (Ts) and do the index/barcode reads shortly after
		runInfo.Run.Reads.Read.each { r ->
			if (!'Y'.equalsIgnoreCase(r.@'IsIndexedRead')) {
				read_structure.add(r.@'NumCycles' + "T")
			}
		}

		computeIndexReadStructure(ssLines, header, "index", index, 1, read_structure, has_umi) // Update the overall read structure with the read structure computed just for the index read.
		computeIndexReadStructure(ssLines, header, "index2", index2, 2, read_structure, has_umi) // Update the overall read structure with the read structure computed just for the index2 read.

		return [num_reads, run_id, flowcell, machine_name, lanecount, read_structure.join(), has_umi]
	}

	static def retrieveListOfTiles(xmlFilePath, instrument, run_directory) {
		def runInfo = new XmlParser().parse(new File(xmlFilePath))
		def tiles = []
		if ('novaseq'.equalsIgnoreCase(instrument)) {
			tiles = runInfo.Run.FlowcellLayout.TileSet.Tiles.Tile
		} else if('miseq'.equalsIgnoreCase(instrument)) {
			tiles = runInfo.Run.TileSelection.Lane.Tile
		} else {
			// Should think about how to handle other instrument types. It could be a NextSeq or Firefly or a new instrument
		}

		def preProcessDir = new File(run_directory + "/Nextflow/preProcess")
		def tilesFile = new File("${preProcessDir}/tilesList.csv")

		if(!preProcessDir.exists()) { preProcessDir.mkdirs()}

		tiles.each {
			r ->
			tilesFile << r.text() + "\n"
		}

		return [tilesFile, preProcessDir]
	}

	static def parseSampleSheet(sampleSheetPath) {
	  // Parse samplesheet.csv and return the experiment name and number of header row to skip
	  def header = 0
	  def experiment = 'NA'
	  def sampleSheetFile = new File(sampleSheetPath)
	  sampleSheetFile.eachLine { line, lineNumber ->
	    def fields = line.split(',')
	    //Randomly check rows for column names to demarcate end of header.
	    if (fields.length>5 && ( fields[0].trim().equalsIgnoreCase('Sample_ID') || fields[1].trim().equalsIgnoreCase('Sample_Name')) ) {
	        header = (lineNumber.toInteger() - 1)
	    }
	    if(fields.length>=2 && line.split(',')[0].trim().equalsIgnoreCase('Experiment Name')) {
	      experiment = line.split(',')[1]
	    }
	  }
	  return [header, experiment]
	}

	static def randomGenerator(n) {
		def alphanumeric_set = (('A'..'Z')+('0'..'9')).join()
		return new Random().with {
	    (1..n).collect { alphanumeric_set[ nextInt( alphanumeric_set.length() ) ] }.join()
	  }
	}

	static def getBarcodeColumns(colnames, targetCol) {
		def targetColId = -1
		colnames.split(",").eachWithIndex { colName, colIdx -> if (colName.trim().equalsIgnoreCase(targetCol)){ targetColId =  colIdx } }
		return targetColId
	}

	static def computeIndexReadStructure(ssData, headerEnd, indexType, colIndex, position, read_structure, has_umi) {
		if(indexType && colIndex >= 0) { // if a column named 'index' exists, this is likely to exist but 'index2' may not always does.
				// NNNNACGT
				def index_sequence = ssData.get(headerEnd.toInteger() + 1).split(",")[colIndex].toUpperCase() // Barcode sequence
				if(index_sequence.startsWith("N") && !index_sequence.endsWith("N") ) {
					has_umi = true // now we know this run has UMI's, update/flip has_umi
					def index_read_structure = (index_sequence.lastIndexOf("N") + 1).toString() + "T"   // Since we know that index sequence has UMI's at the beginnng, the last position of UMI should tell us about the size of the UMI
					index_read_structure += (index_sequence.length() - (index_sequence.lastIndexOf("N") + 1)).toString() + "B"
					read_structure.addAll(position,index_read_structure) // Update the overall read structure with the read structure computed just for the index1 read.

				}
				// ACGTNNNN
				else if(!index_sequence.startsWith("N") && index_sequence.endsWith("N")) {
					has_umi = true
					def index_read_structure = index_sequence.indexOf("N").toString() + "B"
					index_read_structure += (index_sequence.length() - index_sequence.indexOf("N")).toString() + "T"
					read_structure.addAll(position,index_read_structure)

				}
				// NNNNACGTNNNN
				else if(index_sequence.startsWith("N") && index_sequence.endsWith("N")) {
					has_umi = true
					def index_read_structure = index_sequence.findIndexOf{ it !="N"}.toString + "T"
					index_read_structure += ((index_sequence.findLastIndexOf{ it !="N"} - index_sequence.findIndexOf{ it !="N"}) + 1).toString() + "B"
					index_read_structure += (index_sequence.length() - (index_sequence.findIndexOf{ it !="N"} + 1)).toString() + "T"
					read_structure.addAll(position,index_read_structure)

				}
				// No UMI at all AGCTGACT, hence all of the index sequence is a barcode
				else if(index_sequence.findIndexOf{ it =="N"} < 0 ) {
					def index_read_structure = index_sequence.length().toString() + "B"
					read_structure.addAll(position,index_read_structure)
				}
				else { // ACGTNNNNAGCT
					 exit 1, "Nested UMI is not supported!"
				}
		}

	}

	static def retrieveLibraryParams(preProcessDir, samplesheet, lanecount, headerLines, num_reads) {
		// Creates barcodes, multiplex and library params file for use by picard

		for( int laneID in 1 .. lanecount.toInteger()) {

			String laneNo = laneID
			File barcodeParamsFile = new File("${preProcessDir}/lane${laneNo}_barcode_params.txt")
			File multiplexParamsFile = new File("${preProcessDir}/lane${laneNo}_multiplex_params.txt")
			File libraryParamsFile = new File("${preProcessDir}/lane${laneNo}_library_params.txt") // Used only to generate unaligned BAMs, not expensive to create this file regardless
			if (num_reads > 3) {
				barcodeParamsFile << "barcode_sequence_1\tbarcode_sequence_2\tbarcode_name\tlibrary_name\n"
				multiplexParamsFile << "OUTPUT_PREFIX\tBARCODE_1\tBARCODE_2\n"
				multiplexParamsFile << "L${laneNo}_unassigned\tN\tN\n"
				libraryParamsFile << "OUTPUT\tBARCODE_1\tBARCODE_2\tLIBRARY_NAME\tSAMPLE_ALIAS\n"
				libraryParamsFile << "L${laneNo}_unassigned.bam\tN\tN\tunassigned\tunassigned\n"
				Channel.fromPath(samplesheet)
					.splitCsv(header: true, skip: headerLines, strip: true)
					.map { row ->
							barcodeParamsFile << "${row.index}\t${row.index2}\t${row.I7_Index_ID}+${row.I5_Index_ID}\t${MiscUtils.reformatSampleId(row.Sample_ID)}\n"
							multiplexParamsFile << "L${laneNo}_${MiscUtils.reformatSampleId(row.Sample_ID)}\t${row.index}\t${row.index2}\n"
							libraryParamsFile << "L${laneNo}_${MiscUtils.reformatSampleId(row.Sample_ID)}.bam\t${row.index}\t${row.index2}\t${MiscUtils.reformatSampleId(row.Sample_ID)}\t${MiscUtils.reformatSampleId(row.Sample_ID)}\n"
					}
			} else {
					barcodeParamsFile << "barcode_sequence_1\tbarcode_name\tlibrary_name\n"
					multiplexParamsFile << "OUTPUT_PREFIX\tBARCODE_1\n"
					multiplexParamsFile << "L${laneNo}_unassigned\tN\n"
					libraryParamsFile << "OUTPUT\tBARCODE_1\tLIBRARY_NAME\tSAMPLE_ALIAS\n"
					libraryParamsFile << "L${laneNo}_unassigned.bam\tN\tunassigned\tunassigned\n"
					Channel.fromPath(samplesheet)
						.splitCsv(header: true, skip: headerLines, strip: true)
						.map { row ->
								barcodeParamsFile << "${row.index}\t${row.I7_Index_ID}\t${MiscUtils.reformatSampleId(row.Sample_ID)}\n"
								multiplexParamsFile << "L${laneNo}_${MiscUtils.reformatSampleId(row.Sample_ID)}\t${row.index}\n"
								libraryParamsFile << "L${laneNo}_${MiscUtils.reformatSampleId(row.Sample_ID)}.bam\t${row.index}\t${MiscUtils.reformatSampleId(row.Sample_ID)}\t${MiscUtils.reformatSampleId(row.Sample_ID)}\n"
						}
			}
		}
	}

	static def getReadGroups(samplesheet, headerLines, num_reads, lanecount) {
		// Creates a channel comprised of library name and their read group ids
		 
		 def unalignedRdg  = Channel.from([['unassigned', 'unassigned:N:']]) // read group for unaligned BAM, only Barcode 1

		 if(num_reads > 3) { unalignedRdg  = Channel.from([['unassigned', 'unassigned:N:N']]) } // It should not really matter but the format is this when there is a second barcode
		
		 // def x = Channel.from([['hi', 'low']])
		 def readGroups = Channel.fromPath(samplesheet)
			.splitCsv(header: true, skip: headerLines, strip: true)
			.map { row ->
				if (num_reads > 3) {
					[MiscUtils.reformatSampleId(row.Sample_ID), "${MiscUtils.reformatSampleId(row.Sample_ID)}:${row.index}:${row.index2}"]
				} else {
					[MiscUtils.reformatSampleId(row.Sample_ID), "${MiscUtils.reformatSampleId(row.Sample_ID)}:${row.index}:"]
				}
				
			}
			
			readGroups.mix(unalignedRdg)
	}
}
