import java.io.IOException;

import chemaxon.struc.Molecule;
import chemaxon.formats.MolImporter;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.marvin.calculations.pKaPlugin;
import chemaxon.marvin.plugin.PluginException;

import java.text.DecimalFormat;

public final class ChargeMicrospecies_pKa {

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

	pKaPlugin plugin = new pKaPlugin();
	plugin.setpH(7.0);

	Molecule frags[] = OriginalMol.convertToFrags();
	for(int frag_count = 0; frag_count < frags.length; frag_count++){

	    //Retrieve molecule, set implicit hydrogens and aromatize
	    Molecule fragment = frags[frag_count];

	    try{
		plugin.setMolecule(fragment);
		plugin.run();
		
		int msCount = plugin.getMsCount();
		for (int msIndex = 0; msIndex < msCount; msIndex++) {
		    
		    Molecule ms = plugin.getMsMolecule(msIndex);
		    double msDistr = plugin.getSingleMsDistribution(msIndex);
		    
		    // optional step: skip microspecies if its distribution is too low
		    if (msDistr < 0.001) {
			continue;
		    }
		    
		    // write structures in dearomatized form (-a option)
		    // the plugin generates microspecies in aromatized form by default
		    String struct = MolExporter.exportToFormat(ms, "smiles:-a");
		    System.out.println(frag_count + "\t" + msDistr + "\t" + struct);
		    
		}
	    }catch(PluginException PIE){
		    System.out.println("Error: "+PIE.getMessage());
	    }catch(IOException IOE) {
		System.out.println("Error: "+IOE.getMessage());
	    }
	}
    }
}
