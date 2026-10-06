import java.io.IOException;
import java.util.Arrays;
import chemaxon.formats.MolImporter;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.struc.Molecule;
import chemaxon.marvin.calculations.pKaPlugin;
import chemaxon.marvin.calculations.MajorMicrospeciesPlugin;
import chemaxon.marvin.calculations.TautomerizationPlugin;
import chemaxon.marvin.plugin.PluginException;

// Some documentation
// https://docs.chemaxon.com/display/docs/pKa_calculation.html
// Code examples from
// https://apidocs.chemaxon.com/jchem/doc/dev/java/api/chemaxon/marvin/calculations/pKaPlugin.html
public final class pKaValues {

    public static void main(String[] args){

	//Import Mol File
	Molecule mol = new Molecule();
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

	//Ready Plugins
	pKaPlugin pkaplugin = new pKaPlugin();
	MajorMicrospeciesPlugin mmsPlugin = new MajorMicrospeciesPlugin();
	TautomerizationPlugin tPlugin = new TautomerizationPlugin();
	tPlugin.setDominantTautomerDistributionCalculation(true);
	
	try { 
	    pkaplugin.setMolecule(mol);
	    pkaplugin.run();

	    double[] pHValues = pkaplugin.getMacropKaValues(pKaPlugin.ACIDIC);
	    for (int i = 0; i < pHValues.length; ++i) {
		tPlugin.setMolecule(mol);
		tPlugin.setpH(pHValues[i]);
		tPlugin.run();
		Molecule tautomer = tPlugin.getStructure(0);

		mmsPlugin.setMolecule(tautomer);
		mmsPlugin.setpH(pHValues[i]);
		mmsPlugin.run();
		Molecule microspecies = mmsPlugin.getMajorMicrospecies();

		try{
		    String MolString = MolExporter.exportToFormat(microspecies, "inchi:SAbs");
		    System.out.println(String.valueOf(pHValues[i])+"\t"+MolString);
		}catch(IOException IOE) {
		    System.out.println("Error: "+IOE.getMessage());
		}
	    }
	}catch(PluginException PIE){
	    System.out.println("Error: "+PIE.getMessage());
	}
    }
}


		
	    // double[][] MicroSpecies = plugin.getMsDistributions();
	    // System.out.println(String.valueOf(MicroSpecies.length)+"\t"+String.valueOf(Values.length));
	    
		// get microspecies data (molecule and distribution)
		// int count = plugin.getMsCount();
		// for (int j=0; j < count; ++j) {
		//    Molecule ms = plugin.getMsMolecule(j);
		//    double distr = plugin.getSingleMsDistribution(j);
		//    System.out.println(String.valueOf(pHValues[i])+"\t"+String.valueOf(distr));
		// }
		// System.out.println(String.valueOf(Values[i]));
		// System.out.println(String.valueOf(MicroSpecies[i].length));
//	    }
