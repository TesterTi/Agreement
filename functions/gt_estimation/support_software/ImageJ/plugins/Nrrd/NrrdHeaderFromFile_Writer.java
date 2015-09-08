// NrrdHeaderFromFile_Writer
// Plugin to write out a nrrd detached header for a specified
// amiramesh or Biorad PIC file

// (c) Gregory Jefferis 2007
// Department of Zoology, University of Cambridge
// jefferis@gmail.com
// http://flybrain.stanford.edu/nrrd/
// All rights reserved
// Source code released under Lesser Gnu Public License v2

// Use this as an example of how to write detached nrrd headers with other 
// ImageJ writer plugins.  To underline the specific code changes required
// I have decided to extend the Biorad_Writer class for this example.
// However an alternative would have been to edit Biorad_Writer directly
// and provide a method for writing the detached header. 

// Compiling:
// NB this needs Nrrd_Writer, Amira_Reader and Biorad_Reader to be in the
// class path.  Easiest way to do that is to have all the .class files in
// the same directory.  Note also that Nrrd_Writer.java contains the
// definition for NrddFileInfo which must also be in the class path.

import ij.IJ;
import ij.io.FileInfo;
import ij.io.OpenDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;

import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class NrrdHeaderFromFile_Writer implements PlugIn {
	private boolean interactiveMode;
	
	public void run(String arg) {
		// if an argument is supplied, consider that
		// we are in batch mode and suppress messages
		if (!arg.equals("")) {
			interactiveMode=false;
		}
		
		OpenDialog od = new OpenDialog("Select Input Image ...", arg);
		String directory = od.getDirectory();
		String fileName = od.getFileName();
		if (fileName==null)
			return;
		IJ.showStatus("Parsing file information for: " + directory + fileName);		
		

		// Read first 132 bytes to heck file magic
		byte[] buf = new byte[132];
		try {
			InputStream is = new FileInputStream(directory + fileName);
			is.read(buf, 0, 132);
			is.close();
		} catch (IOException e) {
			// Couldn't open the file for reading
			return;
		}
		// Fetch first line of file
		String firstLine="";
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(buf),"US-ASCII"));
			firstLine=br.readLine();
		} catch (Exception e) {
		}
		
		try {
			String name=fileName.toLowerCase();
			Matcher amiraMeshDef = Pattern.compile("#\\s+AmiraMesh.*?(BINARY|ASCII)").matcher(firstLine);			
			if(amiraMeshDef.find() || name.endsWith(".am.gz")) {
				parseAmiraFile(directory, fileName);
				return;
			}
			
			if(name.endsWith(".pic.gz") || ( buf[54]==57 && buf[55]==48) ) {
				parseBioradFile(directory, fileName);
				return;
			}
			if (name.endsWith(".bin") || name.endsWith(".bin.gz") ) {
				parseTorstenRawFile(directory, fileName);
				return;
			}
		}
		catch (Exception e) {
			IJ.showStatus("");
			if(interactiveMode) {
				IJ.showMessage("NrrdHeaderFromFile_Writer:", ""+e);
			} else {
				IJ.log("NrrdHeaderFromFile_Writer: "+e);
			}
			return;
		}
	}
	void writeNrrdHeader(FileInfo fi, Calibration cal,String nrrdFilePath) throws IOException {
		// Figure out if the header and data files are in the same directory
		String nrrdHeader;
		if(fi.directory.equals((new File(nrrdFilePath)).getParent())){
			// Yes same directory
			nrrdHeader=Nrrd_Writer.makeDetachedHeader(fi,cal,true);			
		} else {	
			// No, have to write full file path to data file
			// would be nice to (have the option to) write a relative path
			// but there's no built in way to do that
			nrrdHeader=Nrrd_Writer.makeDetachedHeader(fi,cal,false);
			nrrdHeader=nrrdHeader+"data file: "+new File(fi.directory,fi.fileName).getPath()+"\n";
		}
		FileOutputStream out = new FileOutputStream(nrrdFilePath);
		Writer bw = new BufferedWriter(new OutputStreamWriter(out));
		bw.write(nrrdHeader+"\n");
		bw.close();
		IJ.showStatus("Saved detached nrrd header "+nrrdFilePath);	
	}
	void writeNrrdHeader(FileInfo fi, Calibration cal) throws IOException {
		// Set up the detached nrrd header
		String nrrdFileName=fi.fileName+".nhdr";
		writeNrrdHeader(fi, cal, (new File(fi.directory,nrrdFileName)).getPath());
	}
	
	boolean parseAmiraFile(String directory, String fileName) throws IOException {
		FileInfo fi = null; Calibration cal = null;
		
		Amira_Reader amr = new Amira_Reader();
		try {
			fi = amr.getHeaderInfo(directory, fileName);
			cal= amr.getCalibration();
		} catch (IOException e) {
			return false;
		}		
		writeNrrdHeader(fi,cal);
		return true;
	}
	boolean parseTorstenRawFile(String directory, String fileName) throws IOException {
		FileInfo nfi = null; Calibration cal = null;
		
		TorstenRaw_GZ_Reader trr = new TorstenRaw_GZ_Reader();
		try {
			nfi = trr.getHeaderInfo(directory, fileName);
			cal= trr.getCalibration();
		} catch (IOException e) {
			return false;
		}
		if(nfi.fileName.toLowerCase().endsWith(".gz")) nfi.compression=NrrdFileInfo.GZIP;
		String nrrdFileName=nfi.directory+".nhdr";
		writeNrrdHeader(nfi,cal,nrrdFileName);
		return true;
	}
	boolean parseBioradFile(String directory, String fileName) throws IOException {
		FileInfo fi = null; Calibration cal=null;
		
		GJ_Biorad_Reader br = new GJ_Biorad_Reader();
		try {
			br.setFile(directory,fileName);
			
			fi = br.getHeaderInfo();
			if(fileName.toLowerCase().endsWith(".gz")) fi.compression = NrrdFileInfo.GZIP;
			if (IJ.debugMode) IJ.log("FileInfo: "+fi);
			cal= br.getBioRadCalibration();
		} catch (IOException e) {
			return false;
		}		

		writeNrrdHeader(fi, cal);
		return true;
	}
	
}
