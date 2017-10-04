import java.io.IOException;

import chemaxon.formats.MolImporter;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.struc.Molecule;
import chemaxon.marvin.calculations.MajorMicrospeciesPlugin;
import chemaxon.marvin.calculations.TautomerizationPlugin;
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

	//Ready Charge Plugin
	MajorMicrospeciesPlugin mms_plugin = new MajorMicrospeciesPlugin();
	mms_plugin.setpH(pH);

	//Ready Tautomerizer Plugin
	TautomerizationPlugin t_plugin = new TautomerizationPlugin();
	t_plugin.setDominantTautomerDistributionCalculation(true);
	t_plugin.setpH(pH);

	//Molecule for fusing all fragments
	Molecule FusedMol = new Molecule();

	//Fragment Molecule and iterate
	Molecule frags[] = mol.convertToFrags();
	for(int i = 0; i < frags.length; i++){
	    
	    try{

		//Run Charge Plugin
		mms_plugin.setMolecule(frags[i]);
		mms_plugin.standardize(frags[i]);
		mms_plugin.run();
		Molecule mms_mol = mms_plugin.getMajorMicrospecies();
		
		//Run Tautomer Plugin
		t_plugin.setMolecule(mms_mol);
		t_plugin.standardize(mms_mol);
		t_plugin.run();
		Molecule t_mol = t_plugin.getStructure(0);
		
		//Fuse fragment
		FusedMol.fuse(t_mol,false);
		
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