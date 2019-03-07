DIR ="W:/Analysis/"; //Folder where raw data is 
myDIR = newArray("AJ0160_P1"); //Name of OmeTif w/o ome.tif 
myType =newArray("/CroppedData/", "/FullStacks/");
suffix = ".tif";
for (k=0; k< myDIR.length; k++) {
	for (j=0; j< myType.length; j++) {
		myDIRin = DIR + myDIR[k] + myType[j];
		print(myDIRin);
		list = getFileList(myDIRin);
		for (i=0; i < list.length; i++) {
			if (endsWith(list[i], suffix)) {
		// Make sure to change channel below 
			open(myDIRin+"\\"+list[i]);
			run("Properties...", "channels=32 slices=1 frames=1 unit=pixel pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
			run("Save");	
			run("Close");
			}
		}
	}
}




