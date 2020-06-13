package genepub;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.StringTokenizer;

/**
 * 
 * 
 * @author mahogny
 *
 */
public class ReduceFastaByGenes {
	
	
	public static ArrayList<String> readGeneList() throws IOException{
		System.out.println("Reading gene list");
		ArrayList<String> listGenes=new ArrayList<String>();
		BufferedReader br=new BufferedReader(new FileReader("/home/mahogny/Dropbox/applyPI/umeå/project/bias/genes_18k.csv"));
		String line;
		while((line=br.readLine())!=null) {
			listGenes.add(line);
		}
		br.close();
		
		return listGenes;
	}
	
	public static void main(String[] args) throws NumberFormatException, IOException {
		
		System.out.println("Reading gene list");
		HashSet<String> listGenes=new HashSet<String>(readGeneList());

		
		//Typical line
		//>ENSMUST00000177564.1 cdna chromosome:GRCm38:14:54122226:54122241:1 gene:ENSMUSG00000096176.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:Trdd2 description:T cell receptor delta diversity 2 [Source:MGI Symbol;Acc:MGI:4439546]

		//parse out gene_symbol:Trdd2
		
		
		PrintWriter pw=new PrintWriter(new File("/home/mahogny/Dropbox/applyPI/umeå/project/bias/red_cdna.fa"));
		BufferedReader br=new BufferedReader(new FileReader("/home/mahogny/Dropbox/applyPI/umeå/project/bias/Mus_musculus.GRCm38.cdna.all.fa"));
		String line=br.readLine();
		while(line!=null) {
			String fastaheader=line;
			
			boolean keep=false;
			
			int index=line.indexOf("gene_symbol:");
			if(index!=-1) {
				line=line.substring(index+"gene_symbol:".length());
				index=line.indexOf(' ');
				if(index!=-1) {
					line=line.substring(0,index);
				}
				String geneSymbol=line;
				if(listGenes.contains(geneSymbol)) {
					keep=true;
				}
			}
			
			
			if(keep) {
				pw.println(fastaheader);
				line=br.readLine();
				while(line!=null && !line.startsWith(">")) {
					pw.println(line);
					line=br.readLine();
				}
			} else {
				line=br.readLine();
				while(line!=null && !line.startsWith(">")) {
					line=br.readLine();
				}
			}
		}
		br.close();
		
		
	}

}
