package genepub;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;

public class ReadCOSMIC {

	
	
	public static void main(String[] args) throws IOException {

		
		HashMap<String, Integer> geneCount=new HashMap<String, Integer>();
		HashMap<String, Integer> geneLen=new HashMap<String, Integer>();
		
		BufferedReader br=new BufferedReader(
				new InputStreamReader(
						new GZIPInputStream(
								new FileInputStream(
										new File("/home/mahogny/umeå/project/bias/input/cosmic/CosmicGenomeScreensMutantExport.tsv.gz")))));

		int i=0;
		
		String line;
		br.readLine();
		while((line=br.readLine())!=null) {
			i++;
			
			if(i%1000000==0) {
				System.out.println(i);
				//break;
			}
			//System.out.println(line);
			
			StringTokenizer stok=new StringTokenizer(line,"\t");
			stok.nextToken();
			String enst=stok.nextToken();
			//System.out.println(enst);
			
			enst=enst.substring(0, enst.indexOf("."));
			
			Integer v=geneCount.get(enst);
			if(v==null)
				v=1;
			else
				v++;
			geneCount.put(enst, v);
			
			
			int len=Integer.parseInt(stok.nextToken());
			
			geneLen.put(enst, len);
			
			//System.exit(0);
			
		}
		
		
		br.close();
		
		System.out.println("done reading");
		
		PrintWriter pw=new PrintWriter(new File("/home/mahogny/umeå/project/bias/input/cosmic/genecount.csv"));
		pw.println("gene,count,length");
		for(String g:geneCount.keySet())
			pw.println(g+","+geneCount.get(g)+","+geneLen.get(g));
		pw.close();
		System.out.println("done");
		
	}
}
