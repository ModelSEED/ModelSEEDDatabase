import java.io.IOException;

import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.formats.MolImporter;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.marvin.calculations.MajorMicrospeciesPlugin;
import chemaxon.marvin.calculations.TautomerizationPlugin;
import chemaxon.calculations.hydrogenize.Hydrogenize;
import chemaxon.marvin.plugin.PluginException;

import java.text.DecimalFormat;

public final class AtomNumbers {

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

	Molecule frags[] = OriginalMol.convertToFrags();
	for(int frag_count = 0; frag_count < frags.length; frag_count++){

	    //Retrieve molecule, set implicit hydrogens and aromatize
	    Molecule fragment = frags[frag_count];

	    int atom_total = fragment.getAtomCount();
	    for(int atom_count = 0; atom_count < atom_total; atom_count++){
		MolAtom atom = fragment.getAtom(atom_count);
		System.out.println(frag_count+"\t"+atom_count+"\t"+atom.getSymbol());
	    }
	    try{

		String InChIString = MolExporter.exportToFormat(fragment, "smiles");
		System.out.println("Wotcha");
		
	    }catch(IOException IOE) {
		System.out.println("Error: "+IOE.getMessage());
	    }
	}
    }
}
