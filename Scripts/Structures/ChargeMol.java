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

	// Hydrogenize plugin
	Hydrogenize hydro = new Hydrogenize();

	// MMS Plugin
	MajorMicrospeciesPlugin mmsPlugin = new MajorMicrospeciesPlugin();
	mmsPlugin.setpH(pH);

	// Tautomer plugin
	// Set it so that we retrieve the dominant tautomer at pH of 7
	// Precise calculations may take time so set to five seconds
	TautomerizationPlugin tPlugin = new TautomerizationPlugin();
	tPlugin.setpH(pH);
	tPlugin.setDominantTautomerDistributionCalculation(true);

	// Molecule for fusing all fragments
	Molecule FusedMol = new Molecule();

	// Fragment Molecule and iterate
	// For each fragment we:
	// i)   convert to implicit hydrogens
	// ii)  aromatize
	// iii) find dominant tautomer at pH7
	// iv)  find major microspecies at pH7
	// v)   convert to explicit hydrogens
	Molecule frags[] = OriginalMol.convertToFrags();
	for(int i = 0; i < frags.length; i++){

		//Retrieve molecule, set implicit hydrogens and aromatize
		Molecule fragment = frags[i];
		hydro.convertExplicitHToImplicit(fragment);
		fragment.aromatize();

		try{
		    
		    // Run Tautomer Plugin
		    tPlugin.setMolecule(fragment);
		    tPlugin.run();

		    // Dominant tautomer is the first one
		    Molecule tautomer = tPlugin.getStructure(0);

		    // Run MMS Plugin
		    mmsPlugin.setMolecule(tautomer);
		    mmsPlugin.run();

		    //Get Major microspecies as molecule
		    Molecule microspecies = mmsPlugin.getMajorMicrospecies();

		    // Convert to explicit hydrogens
		    hydro.convertImplicitHToExplicit(microspecies);
		    
		    //Fuse fragment
		    FusedMol.fuse(microspecies,false);
		
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
