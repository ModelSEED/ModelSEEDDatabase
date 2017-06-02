import java.io.IOException;

import chemaxon.formats.MolImporter;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.struc.Molecule;
import chemaxon.marvin.calculations.MajorMicrospeciesPlugin;
import chemaxon.marvin.plugin.PluginException;

public final class ChargeMol {

    public static void main(String[] args){
	Molecule mol = new Molecule();
	double pH = 7.0;

	//Import Mol File
	try {
	
	    MolImporter mi = new MolImporter(args[0]);

	    try {
		mol = mi.read();
	    } finally {
		mi.close();
	    }

	} catch (MolFormatException e) {
	    throw new IllegalArgumentException("Invalid molecule format", e);
	} catch (IOException e) {
	    throw new IllegalArgumentException("Error reading input file", e);
	}

	//Check for pH
	if(args.length == 2){
	    pH = Double.parseDouble(args[1]);
	}

	//Ready Plugin
	MajorMicrospeciesPlugin plugin = new MajorMicrospeciesPlugin();
	//System.out.println("pH "+String.valueOf(pH));
	plugin.setpH(pH);

	//Molecule for fusing all fragments
	Molecule FusedMol = new Molecule();

	//Fragment Molecule and iterate
	Molecule frags[] = mol.convertToFrags();
	for(int i = 0; i < frags.length; i++){
	    
	    try{
		
		//Run Pluging
		plugin.setMolecule(frags[i]);
		plugin.run();
		
		//Get Major microspecies as molecule
		Molecule FragMMS = plugin.getMajorMicrospecies();
		
		//Fuse fragment
		FusedMol.fuse(FragMMS,false);
		
	    }catch(PluginException PIE){
		System.out.println("Error: "+PIE.getMessage());
	    }

	    //try{
	
		//String FragString = MolExporter.exportToFormat(frags[i], "inchi");
		//String FragMMSString = MolExporter.exportToFormat(FragMMS, "inchi");
		//System.out.println("Fragment "+Integer.toString(i)+" : "+FragString);
		//System.out.println("Fragment "+Integer.toString(i)+" : "+FragMMSString);
		
	    //}catch(IOException IOE) {
	    //System.out.println("Error: "+IOE.getMessage());
	    //}
	}
	
	//Dearomatize
	try{
	    
	    FusedMol.dearomatize();

	}catch(SecurityException SE){
	    System.out.println("Error: "+SE.getMessage());
	}

	//Print to string
	try{
	
	    String FusedString = MolExporter.exportToFormat(FusedMol, "mol");
	    System.out.print(FusedString);

	}catch(IOException IOE) {
	    System.out.println("Error: "+IOE.getMessage());
	}
    }
}