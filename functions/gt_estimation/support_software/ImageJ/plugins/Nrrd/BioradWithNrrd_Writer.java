// BioradWithNrrd_Writer
// ---------------------
// Plugin to write out an ImagePlus as Biorad PIC 
// and at the same time write out a nrrd detached header
// See http://teem.sourceforge.net/nrrd/
// and http://flybrain.stanford.edu/nrrd/

// (c) Gregory Jefferis 2007
// Department of Zoology, University of Cambridge
// jefferis@gmail.com
// All rights reserved
// Source code released under Lesser Gnu Public License v2

// Use this as an example of how to write detached nrrd headers with other 
// ImageJ writer plugins.  To underline the specific code changes required
// I have decided to extend the Biorad_Writer class for this example.
// However an alternative would have been to edit Biorad_Writer directly
// and provide a method for writing the detached header. 

// Compiling:
// NB this needs Nrrd_Writer, Nrrd_Reader and Biorad_Writer to be in the
// class path.  Easiest way to do that is to have all the .class files in
// the same directory.  Note also that Nrrd_Reader.java contains the
// definition for NrddFileInfo which is also required.

import ij.IJ;
import ij.io.FileInfo;
import ij.measure.Calibration;
import ij.plugin.PlugIn;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;

public class BioradWithNrrd_Writer extends Biorad_Writer implements PlugIn {

	void writeImage(FileInfo fi, Calibration cal) throws IOException {
		// Write the Biorad file
		super.writeImage(fi, cal);
		
		// Set data offset
		// NB this is easy for Biorad because offset is fixed
		// in other cases, may be necessary to extract the offset from the base
		// class (eg Amira_Writer), which may complicate things
		// ie may actually be nessessary to modify the base class itself
		fi.offset=76;
		// Set up the detached nrrd header
		String nrrdHeader=Nrrd_Writer.makeDetachedHeader(fi,cal,true);
		String nrrdFileName=fi.fileName+".nhdr";
		
		// Write it out
		FileOutputStream out = new FileOutputStream(new File(fi.directory,nrrdFileName));
		Writer bw = new BufferedWriter(new OutputStreamWriter(out));
		bw.write(nrrdHeader+"\n");
		bw.close();
		IJ.showStatus("Saved detached nrrd header "+nrrdFileName);
	}
}
