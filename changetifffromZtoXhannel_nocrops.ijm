DIR ="W:/Analysis/"; 
myDIR = newArray("AJ0160_P1","AJ0160_P2", "AJ0160_P3", "AJ0160_P5","AJ0176_P1","AJ0176_P2","AJ0176_P3","AJ0176_P4","AJ0176_P5"); 
myType =newArray("/CroppedData/", "/SplitData/");
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




