package genepub;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintWriter;


import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.input.SAXBuilder;
import org.jdom2.input.sax.XMLReaders;
import org.jdom2.output.XMLOutputter;

/**
 * Extract PubMed data as a table. from ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline
 * 
 * @author Johan Henriksson
 *
 */
public class ExtractDate {

	public static void printXML(Element e) {
		System.out.println(new XMLOutputter().outputString(e));
	}

	
	public static void doFile(File fXmlFile) throws IOException {
		SAXBuilder builder = new SAXBuilder(XMLReaders.NONVALIDATING);
        builder.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd", false);

		
		System.out.println("Start "+fXmlFile );

		File outFile=new File(fXmlFile.getParentFile(), fXmlFile.getName()+".csv");
		
		System.out.println(outFile);

		if(!outFile.exists()) {

			try {
				PrintWriter pw=new PrintWriter(outFile);

				/*
				
				DocumentBuilderFactory factory = DocumentBuilderFactory.newDefaultInstance();
				factory.setNamespaceAware( true );
			//	builder = factory.newDocumentBuilder();
				builder.setEntityResolver( new BlankingResolver() );
				Document = builder.build( (new FileInputStream(fXmlFile) );
				*/
						
				/*
			        String baseString = JDom2Util.getInstance().readFile(base.getPath(), charSet);
			        Document baseDoc =
			            builder.build(new InputSource(new BufferedReader(new StringReader(baseString))));
			        */
				
//				Document document = (Document) builder.build(new GZIPInputStream(new FileInputStream(fXmlFile)));
				Document document = (Document) builder.build((new FileInputStream(fXmlFile)));
				Element eRoot = document.getRootElement();

				for (Element eDoc : eRoot.getChildren()) {

					// printXML(eDoc);
					Element ePubmed = eDoc.getChild("PubmedData");

					//Figure out the pubmed ID. think this is the way
					int id = -1;
					Element eIdList = ePubmed.getChild("ArticleIdList");
					for (Element eId : eIdList.getChildren()) {
						if (eId.getAttributeValue("IdType").equals("pubmed")) {
							id = Integer.parseInt(eId.getText());
						}
					}

					Element eHistory = ePubmed.getChild("History");
					for (Element eDate : eHistory.getChildren()) {
						// Loop over all the PubMedPubDate elements. Take the first one, whatever it is

						double year = Double.parseDouble(eDate.getChild("Year").getText());
						double month = Double.parseDouble(eDate.getChild("Month").getText());
						double day = Double.parseDouble(eDate.getChild("Day").getText());

						double comb = year + month / 12.0 + day / 365.0;

						pw.println(id + "\t" + year + "\t" + month + "\t" + day + "\t" + comb);
						
						//Only take the first ID
						break; 
					}

				}
				pw.close();

			} catch (Exception e) {
				e.printStackTrace();
			}

		}
		
	}
	
	public static void main(String[] args) throws IOException {

//		if(!f.getName().startsWith(".") && f.getName().endsWith(".gz")) {

		//File fXmlFile = new File("pubmed19n0001.xml.gz");
		for(File f:new File("D:\\Henlab\\pubmed_baseline").listFiles()) {
			if(!f.getName().startsWith(".") && !f.getName().endsWith(".csv")) {
				doFile(f);
			}
		}
		System.out.println("done all files");
	}

}

