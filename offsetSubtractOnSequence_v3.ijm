// please download the latest thunderSTORM https://github.com/zitmen/thunderstorm/releases, (solved Issue with calls from macro #19, https://github.com/zitmen/thunderstorm/issues/19)
// 170108 TD v2 save all offset subtracted folder in a separated folder
// 180622 TD v2_1 sequential data sets
// 190624 TD offsetSubtractOnSequence_v2 added offset calculation part and some comments for parameters

macro "Run offset subtraction" { 
	print("----start 'Run offset subtraction macro'----" );

	///// Parameter start

	// File path of your folder which includes your raw subfolders. Use '\\' as file separator
	// folderPath = "G:\\RawData\\190729_2158BBHPTara laser irradiation tests"
	// folderPath = "E:\\DocuBak\\TD_WZ TAB RoSEO\\Tara Porter\\190820_2158BBP increased damage"
	folderPath = "K:\\WZ 1 data\\210902_2158BBP LL_Nile Red_ultra pure H2O (NR in 1X PBS)"

	// Prefix of data name
	dataPrefix = "Data_"
	darkPrefix = "Dark_"

	// data index of offset image sequence
	offsetDataInd = "1";

	// the first data index you want to offset subtract
	dataStart = 9;
	// the final data index you want to offset subtract
	dataEnd = 9;

	// file format
	fmt = ".tif";

	///// Parameter end

	offsetDataName = darkPrefix + offsetDataInd;
	offsetDataFolder = folderPath + "\\" + offsetDataName + "\\";
	offsetSaveName = "AVG_" + offsetDataName;
	// open image sequence
    run("Image Sequence...", "open=[offsetDataFolder] sort");
    run("Z Project...", "projection=[Average Intensity]");
    saveAs("Tiff", folderPath + "\\" + offsetSaveName + fmt);
    while (nImages>0){
			selectImage(nImages);
			close();
    }
	
	// create a separated folder for saving offset subtracted images
	newFolder = folderPath + "\\offsetSubtracted";
	File.makeDirectory(newFolder)

	//for(i=0; i<dataIndex.length; i++) {
    for(i=dataStart; i<=dataEnd; i++) {
		dataIndex = i;				
		
		// open offset image
		open(folderPath + "\\" + offsetSaveName + fmt);
		
		dataName = dataPrefix + dataIndex;
	    print(dataName);
	    //runMacro("C:\\Users\\tding\\My Works\\Project\\5. Transient amyloid binding imaging\\Analysis\\ImageJgetTime.ijm");
	    
	    // create a new folder for saving offset subtracted images
	    newSubFolder = newFolder + "\\" + dataName;
	    File.makeDirectory(newSubFolder)
	    newName = newSubFolder + "\\" + dataName + "0000" + fmt;
	    
	    dataFolder = folderPath + "\\" + dataName + "\\";
	    // open image sequence
    	run("Image Sequence...", "open=[dataFolder] sort");

    	// offset subtraction
		imageCalculator("Subtract create stack", dataName, offsetSaveName+fmt);
		selectWindow("Result of " + dataName);
		
		// save offset subtracted images as an image sequence
		run("Image Sequence... ", "format=TIFF name=[" + dataName + "] save=[" + newName + "]");

		while (nImages>0){
			selectImage(nImages);
			close();
		}			
		print("------");	
	}
	print("----end of 'Run offset subtraction macro'----" );
	//runMacro("C:\\Users\\tding\\My Works\\Project\\5. Transient amyloid binding imaging\\Analysis\\ImageJgetTime.ijm");
}
