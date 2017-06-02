import java.io.IOException;
//import java.util.ArrayList;
//import java.util.List;

import chemaxon.formats.MolImporter;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.struc.Molecule;
//import system.out;

public final class ReadMolFile {

    public static Molecule importMol(String fullPath) {
	try {
	
	    MolImporter mi = new MolImporter(fullPath);

	    try {
		return mi.read();
	    } finally {
		mi.close();
	    }
	} catch (MolFormatException e) {
	    throw new IllegalArgumentException("Invalid molecule format", e);
	} catch (IOException e) {
	    throw new IllegalArgumentException("Error reading input file", e);
	}
    }

    public static void main(String[] args){
	Molecule mol = importMol(args[0]);
	Molecule frags[] = mol.convertToFrags();

	for(int i = 0; i < frags.length; i++){
	    try {
		String smiles = MolExporter.exportToFormat(frags[i], "inchi");
		System.out.println("Fragment "+Integer.toString(i)+" : "+smiles);
		
	    }catch (IOException e) {
		throw new IllegalArgumentException("Error", e);
	    }
	}
    }
}