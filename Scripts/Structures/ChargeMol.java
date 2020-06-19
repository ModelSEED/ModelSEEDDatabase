import java.io.IOException;

import chemaxon.struc.Molecule;
import chemaxon.formats.MolImporter;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.marvin.calculations.MajorMicrospeciesPlugin;
import chemaxon.marvin.calculations.TautomerizationPlugin;
import chemaxon.calculations.hydrogenize.Hydrogenize;
import chemaxon.marvin.plugin.PluginException;

public final class ChargeMol {

    public static void main(String[] args){
	Molecule OriginalMol = new Molecule();
	double pH = 7.0;

	//Import Mol File
	try {
	
	    MolImporter mi = new MolImporter(args[0]);

	    try {
		OriginalMol = mi.read();
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

	//Property keys should be restored, work in progress
	//Enumeration keys = myMol.getPropertyKeys();
	//while (keys.hasMoreElements()) {
	//    String key = (String) keys.nextElement();
	//    String value = myMol.getProperty(key);
	//    for (int i=0; i < molArray.length; ++i) {     
	//	molArray[i].setProperty(key, value);
	//    }
	//}

	//Ready Plugin
	MajorMicrospeciesPlugin mmsPlugin = new MajorMicrospeciesPlugin();
	mmsPlugin.setpH(pH);

	//Ready Tautomer plugin
	//Set it so that we retrieve the dominant normal canonical tautomer at pH of 7
	//Precise calculations may take time so set to five seconds
	TautomerizationPlugin tPlugin = new TautomerizationPlugin();
	tPlugin.setTakeCanonicalForm(true);

	Hydrogenize hydro = new Hydrogenize();

	//Molecule for fusing all fragments
	Molecule FusedMol = new Molecule();

	//Fragment Molecule and iterate
	Molecule frags[] = OriginalMol.convertToFrags();
	for(int i = 0; i < frags.length; i++){
	    
	    try{

		//Retrieve molecule, set implicit hydrogens and aromatize
		Molecule Frag = frags[i];
		hydro.convertExplicitHToImplicit(Frag);
		Frag.aromatize();

		//Run MMS Plugin
		mmsPlugin.standardize(Frag);
		mmsPlugin.setMolecule(Frag);
		mmsPlugin.run();

		//Get Major microspecies as molecule
		Molecule FragMMS = mmsPlugin.getMajorMicrospecies();

		//Run Tautomer Plugin
		tPlugin.standardize(FragMMS);
		tPlugin.setMolecule(FragMMS);
		tPlugin.run();

		//Dominant tautomer is the first one
		Molecule FragMMSTaut = tPlugin.getStructure(0);
		
		//Dearomatize
		FragMMSTaut.dearomatize();

		//Fuse fragment
		FusedMol.fuse(FragMMSTaut,false);
		
	    }catch(PluginException PIE){
		System.out.println("Error: "+PIE.getMessage());
	    }
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