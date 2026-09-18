import java.io.IOException;

import chemaxon.struc.Molecule;
import chemaxon.formats.MolImporter;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.marvin.calculations.MajorMicrospeciesPlugin;
import chemaxon.marvin.calculations.TautomerizationPlugin;
import chemaxon.calculations.hydrogenize.Hydrogenize;
import chemaxon.marvin.plugin.PluginException;

import java.text.DecimalFormat;

public final class ChargeMicrospecies {

    private static final DecimalFormat df = new DecimalFormat("0.000000");

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
	// Molecule FusedMol = new Molecule();

	// Fragment Molecule and iterate
	// For each fragment we:
	// i)   convert to implicit hydrogens
	// ii)  aromatize
	// iii) find dominant tautomer at pH7
	// iv)  find major microspecies at pH7
	// v)   convert to explicit hydrogens
	Molecule frags[] = OriginalMol.convertToFrags();
	for(int frag_count = 0; frag_count < frags.length; frag_count++){

		//Retrieve molecule, set implicit hydrogens and aromatize
		Molecule fragment = frags[frag_count];
		// hydro.convertExplicitHToImplicit(fragment);
		// fragment.aromatize();

		try{
		    
		    // Run Tautomer Plugin
		    // Dominant tautomer is the first one
		    // tPlugin.setMolecule(fragment);
		    // tPlugin.run();
		    // Molecule tautomer = tPlugin.getStructure(0);

		    // Run MMS Plugin
		    mmsPlugin.setMolecule(fragment);
		    mmsPlugin.run();

		    //Get Major microspecies as molecule
		    //Molecule microspecies = mmsPlugin.getMajorMicrospecies();

		    int ms_total = mmsPlugin.getMicrospeciesCount();
		    for(int ms_count = 0; ms_count < ms_total; ms_count++){
			Molecule microspecies = mmsPlugin.getSortedMicrospecies(ms_count);
			double distribution = mmsPlugin.getSortedMsDistribution(ms_count);
		    
			// Convert to explicit hydrogens
			// hydro.convertImplicitHToExplicit(microspecies);
		    
			//Fuse fragment
			//FusedMol.fuse(microspecies,false);

			// String InChIString = MolExporter.exportToFormat(microspecies, "inchi:SAbs,AuxNone,Woff");
			String InChIString = MolExporter.exportToFormat(microspecies, "smiles");
			// System.out.println(frag_count+"\t"+ms_count+"\t"+df.format(distribution)+"\t"+InChIString);
			System.out.println(frag_count+"\t"+ms_count+"\t"+distribution+"\t"+InChIString);
		    }
		    
		}catch(PluginException PIE){
		    System.out.println("Error: "+PIE.getMessage());
		}catch(IOException IOE) {
		    System.out.println("Error: "+IOE.getMessage());
		}
	}
    }
}
