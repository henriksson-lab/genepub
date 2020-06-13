package genepub;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

public class LinkHICtads {

	public static String parseQuote(String s) {
		if(s.startsWith("\"")) {
			return s.substring(1, s.length()-1);
		} else
			return s;
	}

	//Chrom -> Pos -> Gene
	public static HashMap<String, TreeMap<Integer, String>> mapPosGene=new HashMap<String, TreeMap<Integer,String>>();

	
	public static ArrayList<String> getGenes(String chrom, int from, int to) {
		
		if(chrom.startsWith("chr")) {
			chrom=chrom.substring("chr".length());
		}
		
		if(chrom.equals("23"))
			return new ArrayList<String>();
		if(chrom.equals("24"))
			return new ArrayList<String>();
		
		if(!mapPosGene.containsKey(chrom)) {
			System.out.println(chrom);
		}
		
		ArrayList<String> genes=new ArrayList<String>(mapPosGene.get(chrom).tailMap(from).headMap(to).values());
		return genes;
	}
	
	public static HashMap<String, Double> mapGenePmid=new HashMap<String, Double>();
			
	
	public static void main(String[] args) throws IOException {

				
				
		//HashMap<String, Integer> geneCount=new HashMap<String, Integer>();
		//HashMap<String, Integer> geneLen=new HashMap<String, Integer>();
		
		System.out.println("Reading coords");
		//"gene","chrom","pos","ensid","genever","rank_pmid"
		//"9930111J21Rik1","11",48962774,"ENSMUSG00000069893",10,1
		BufferedReader br=new BufferedReader(new FileReader(new File("/home/mahogny/umeå/project/bias/plots/data_chromloc.csv")));
		br.readLine();
		String line;
		while((line=br.readLine())!=null) {
			
			StringTokenizer stok=new StringTokenizer(line,",");
			String gene=parseQuote(stok.nextToken());
			String chrom=parseQuote(stok.nextToken());
			String pos=parseQuote(stok.nextToken());
//String ensid=
			parseQuote(stok.nextToken());
//String genever=
			parseQuote(stok.nextToken());
			Double rank_pmid=Double.parseDouble(parseQuote(stok.nextToken()));
			
			TreeMap<Integer, String> chromDat=mapPosGene.get(chrom);
			if(chromDat==null) {
				mapPosGene.put(chrom, chromDat=new TreeMap<Integer, String>());
			}
			
			chromDat.put((int)Double.parseDouble(pos),gene);
			
			mapGenePmid.put(gene, rank_pmid);
			

		}
		br.close();
		
		System.out.println("Read HiC map");
		br=new BufferedReader(new FileReader(new File("/home/mahogny/umeå/project/bias/input/hic/mm10/mESC.Bonev_2017-raw.domains")));


		
		
		
		int curIndex=0;
		PrintWriter pw=new PrintWriter(new File("/home/mahogny/umeå/project/bias/input/hic/mm10/mESC.Bonev_2017-raw.domains_genepairs.csv"));		
		pw.println("gene\ttad_group\trank_tad");
		while((line=br.readLine())!=null) {
			
			StringTokenizer stok=new StringTokenizer(line,"\t");
			String chrom1=stok.nextToken();
			int from1 = Integer.parseInt(stok.nextToken());
			int to1 = Integer.parseInt(stok.nextToken());

			
			ArrayList<String> genes1=getGenes(chrom1, from1, to1);
			
			if(genes1.size()>1) {
//				System.out.println(genes1);
//				System.out.println();
				
				for(String g:genes1) {
					
					int numOther=genes1.size()-1;
					double sum=0;
					for(String otherGene:genes1) {
						if(!otherGene.equals(g))
							sum += mapGenePmid.get(otherGene);
					}
					sum/=numOther;
					
	//				mapGenePmid.
					
					
					pw.println(g+"\tgroup"+curIndex+"\t"+sum);
				}

			}
			
			curIndex++;
		}
		pw.close();
		
		

		
		
/*
		
		
		
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
		
		
		
		System.out.println("done reading");
		
		PrintWriter pw=new PrintWriter(new File("/home/mahogny/umeå/project/bias/input/cosmic/genecount.csv"));
		pw.println("gene,count,length");
		for(String g:geneCount.keySet())
			pw.println(g+","+geneCount.get(g)+","+geneLen.get(g));
		pw.close();
		System.out.println("done");
		*/
		
	}
	
}
