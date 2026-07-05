import java.io.IOException;

import chemaxon.struc.Molecule;
import chemaxon.formats.MolImporter;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;

import chemaxon.marvin.calculations.TautomerizationPlugin;
import chemaxon.marvin.plugin.PluginException;

public final class ListTautomers {

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

	//Ready Tautomer plugin
	//Set it so that we retrieve the dominant normal canonical tautomer at pH of 7
	//Precise calculations may take time so set to five seconds
	TautomerizationPlugin tPlugin = new TautomerizationPlugin();
	tPlugin.setpH(pH);
	tPlugin.setDominantTautomerDistributionCalculation(true);

	//Fragment Molecule and iterate
	Molecule frags[] = OriginalMol.convertToFrags();
	for(int i = 0; i < frags.length; i++){
	    
	    try{
		
		//Run Tautomer Plugin
		tPlugin.setMolecule(frags[i]);
		tPlugin.run();

		System.out.println("---------Tautomers---------");
		int count = tPlugin.getStructureCount();
		for (int j=0; j < count; ++j) {
		    Molecule tautomer = tPlugin.getStructure(j);
		    double distribution = tPlugin.getDominantTautomerDistribution(j);

		    try{
			
			String InChIKey= MolExporter.exportToFormat(tautomer, "inchikey:SAbs");
			System.out.println(InChIKey+"\t"+distribution);

		    }catch(IOException IOE) {
			System.out.println("Error: "+IOE.getMessage());
		    }
		}
		
	    }catch(PluginException PIE){
		System.out.println("Error: "+PIE.getMessage());
	    }
	}
    }
}
