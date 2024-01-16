import java.io.IOException;
import java.util.Arrays;
import chemaxon.formats.MolImporter;
import chemaxon.formats.MolFormatException;
import chemaxon.struc.Molecule;
import chemaxon.marvin.calculations.pKaPlugin;
import chemaxon.marvin.plugin.PluginException;

// Some documentation
// https://docs.chemaxon.com/display/docs/pKa_calculation.html
// Code examples from
// https://apidocs.chemaxon.com/jchem/doc/dev/java/api/chemaxon/marvin/calculations/pKaPlugin.html
public final class pKaMol {

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

	//Ready Plugin
	pKaPlugin plugin = new pKaPlugin();
	double pH = 7.0;
	plugin.setpH(pH);
	
	//Fragment Molecule and iterate
	Molecule frags[] = mol.convertToFrags();
	for(int frag = 0; frag < frags.length; frag++){

	    try { 
		plugin.setMolecule(frags[frag]);
		plugin.run();
		
		// get pKa values for each atom
		int atom_count = frags[frag].getAtomCount();
		for (int atom=0; atom < atom_count; atom++) {
		    
		    // get ACIDIC and BASIC pKa values
		    double [] apka = plugin.getpKaValues(atom, pKaPlugin.ACIDIC);
		    double [] bpka = plugin.getpKaValues(atom, pKaPlugin.BASIC);
		    
		    String apka_str = "null";
		    if (apka != null) {
			apka_str = String.valueOf(apka[0]);
		    }
		    String bpka_str = "null";
		    if (bpka != null) {
			bpka_str = String.valueOf(bpka[0]);
		    }
		    if(apka_str != "null" || bpka_str != "null"){
			String label = frags[frag].getAtom(atom).getSymbol();
			
			String Output = frag+"\t"+(atom+1)+"\t"+label+"\t"+apka_str+"\t"+bpka_str;

			System.out.println(Output);
		    }
		}
	    }catch(PluginException PIE){
		System.out.println("Error: "+PIE.getMessage());
	    }
	}
    }
}
