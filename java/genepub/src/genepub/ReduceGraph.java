package genepub;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;


/**
 * Takes a graph of gene homologies and removes edges by a type of "triangle inequality"
 * 
 * before we removed Gm genes: 
 * # edges before reduction: 3664858
	# edges after reduction: 143170
	# genes involved: 26204
 * 
 * 
 * # edges before reduction: 2957174
 * # edges after reduction: 794286
 * # genes involved: 10279
 * 
 * with protein blastp
 * # genes involved: 17284
 * # edges before reduction: 4253032
# edges after reduction: 2340616
 * 
 * 
 * @author Johan Henriksson
 *
 */
public class ReduceGraph {

	
	public static HashMap<String, HashMap<String, Double>> mapFromTo=new HashMap<String, HashMap<String,Double>>();
	public static HashSet<Edge> listEdge=new HashSet<Edge>(10000000);

	public static void addEdge(String geneFrom, String geneTo, double pident) {
		
		addDir(geneFrom, geneTo,   pident);
		addDir(geneTo,   geneFrom, pident);
		
	}
	

	public static void addDir(String geneFrom, String geneTo, double pident) {
		
		HashMap<String, Double> mapTo=mapFromTo.get(geneFrom);
		if(mapTo==null) {
			mapTo=new HashMap<String, Double>();
			mapFromTo.put(geneFrom, mapTo);
		}
		mapTo.put(geneTo, pident);
	}
	
	public static void buildGraph() {
		mapFromTo.clear();
		for(Edge e:listEdge) {
			addEdge(e.from, e.to, e.pident);
		}
		System.out.println("# genes involved: "+mapFromTo.size());
		System.out.println("# edge "+listEdge.size());
	}

	
	public static void removeFromGraph(Collection<Edge> removeEdge) {
		
		for(Edge e:removeEdge) {
			mapFromTo.get(e.from).remove(e.to);
			mapFromTo.get(e.to).remove(e.from);
		}
		
	}
	
	

	private static void keepKbest() {
		
		listEdge.clear();
		int k=3;
		for(String from:mapFromTo.keySet()) {
			HashMap<String, Double> mapFrom=mapFromTo.get(from);
				
			//Sort edges
			ArrayList<String> listTo=new ArrayList<String>(mapFrom.keySet());
			Collections.sort(listTo, new Comparator<String>() {
				public int compare(String a, String b) {
					return mapFrom.get(a).compareTo(mapFrom.get(b));
				}
			});
			
			//Put edge on to-keep list
			for(int i=0;i<k && i<listTo.size();i++) {
				String to=listTo.get(i);
				listEdge.add(new Edge(from, to, mapFrom.get(to)));
			}		
		}
		System.out.println("Rebuilding graph after k-best");
		buildGraph();
		
	}
	
	public static void buildEdgeListFromGraph() {
		HashSet<Edge> setEdge=new HashSet<Edge>();
		
		for(String node1:mapFromTo.keySet()) {
			HashMap<String, Double> mapFrom1=mapFromTo.get(node1);
			for(String node2:mapFrom1.keySet()) {
				setEdge.add(new Edge(node1, node2, mapFrom1.get(node2)));
			}
		}

		listEdge.clear();
		listEdge.addAll(setEdge);
	}
	
	
	
	
	
	public static void main(String[] args) throws IOException {

		/////////// Read all the edges. Assume a->b and b->a exist
		readEdgeList(new File("/home/mahogny/Dropbox/applyPI/umeå/project/bias/homology_graph.csv"));
		
		/////////// Build graph
		System.out.println("Building graph");
		buildGraph();

				
		/////////// Only keep the top K neighbors for each gene
		System.out.println("Keeping k best");
		keepKbest();
		
		/////////// The triangle inequality reduction...
		doTriangleIneq();


		System.out.println("# genes involved: "+mapFromTo.size());
		System.out.println("# edge "+listEdge.size());

		
		System.out.println("Writing edges to file");
		PrintWriter pw=new PrintWriter(new File("/home/mahogny/Dropbox/applyPI/umeå/project/bias/homology_graph.red.csv"));
		pw.println("from,to,pident");
		for(Edge e:listEdge) {
			if(e.from.compareTo(e.to)<0)
				pw.println(e.from+","+e.to+","+e.pident);
		}
		pw.close();
		System.out.println("Done");
		
	}


	private static void readEdgeList(File file) throws IOException {
		
		System.out.println("Reading edge list");
		BufferedReader br=new BufferedReader(new FileReader(file));
		String line;
		while((line=br.readLine())!=null) {
			
			StringTokenizer stok=new StringTokenizer(line,"\t");
			
			String geneFrom=stok.nextToken();
			String geneTo=stok.nextToken();
			String pident=stok.nextToken();
			
			Edge e1=new Edge(geneFrom, geneTo, Double.parseDouble(pident));
			listEdge.add(e1);

		}
		br.close();
		
	}


	private static void doTriangleIneq() {
		
		System.out.println("Reducing graph");
		ArrayList<Edge> removeEdges=new ArrayList<Edge>(listEdge.size());
		int i=0;
		for(String node1:mapFromTo.keySet()) {
			if(i%1000==0) {
				System.out.println(i);
			}
			i++;
			for(String node2:mapFromTo.get(node1).keySet()) {
				
				HashMap<String, Double> mapFrom1=mapFromTo.get(node1);
				HashMap<String, Double> mapFrom2=mapFromTo.get(node2);

				//What does 1-2 share in terms of third nodes?
				HashSet<String> commonThird=new HashSet<String>(mapFrom1.keySet());
				commonThird.retainAll(mapFrom2.keySet());
				
				//Now need to compare from-to-third. Try to remove the edge from-to  by comparing from-third and to-third
				double weight12=mapFrom1.get(node2);
				loop3: for(String node3:commonThird) {
					
					//Compare the weights
					double weight13=mapFrom1.get(node3);
					double weight23=mapFrom2.get(node3);
					if(
							weight13 > weight12 & 
							weight23 > weight12) {
						removeEdges.add(new Edge(node1, node2, weight12));
						break loop3;
					}
				}
			}
		}

		System.out.println("# edges before reduction: "+listEdge.size());

		removeFromGraph(removeEdges);
		buildEdgeListFromGraph();
		
		System.out.println("# edges after reduction: "+listEdge.size());
		
	}


}
