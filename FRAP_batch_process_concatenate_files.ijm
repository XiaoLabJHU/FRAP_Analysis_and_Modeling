//Hopefully a working program to go through a list of cell files and concatenate pre and post bleach images together

//Ask user to choose the input directory
directory = getDirectory("Choose input directory");
fileList = getFileList(directory);

for (i=0; i<fileList.length; i++) {
	newfolder = directory + fileList[i];
	newFileList = getFileList(newfolder);
	if (newFileList.length>0){
		for (k=0; k<newFileList.length; k++) {
			file = newFileList[k];
			test_prebleach = indexOf(file,"prebleach");
			test_postbleach = indexOf(file,"postbleach");
			if (test_prebleach > 0) {
				prebleach_filename = newFileList[k];
				open(newfolder + file);
			}
			if (test_postbleach > 0) {
				postbleach_filename = newFileList[k];
				open(newfolder + file);
			}
				
		}
		//check to see if "FRAP_stack.tif" exists; if it does- delete it (ImageJ doesn't overwrite apparently)
		if (File.exists(newfolder + "FRAP_stack.tif") == 1) {
			File.delete(newfolder + "FRAP_stack.tif")
		}
		//concatenate the open images together
		selectWindow(postbleach_filename);
		selectWindow(prebleach_filename);
		run("Concatenate...", "  title=[FRAP_stack] image1="+prebleach_filename+" image2="+postbleach_filename+" image3=[-- None --]");
		saveAs("Tiff", newfolder + "FRAP_stack.tif");
		close();
	}
}
