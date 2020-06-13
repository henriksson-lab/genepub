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

//import javax.xml.parsers.DocumentBuilderFactory;


 
import org.jdom2.Attribute;
import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.JDOMException;
import org.jdom2.filter.Filters;
import org.jdom2.input.DOMBuilder;
import org.jdom2.input.SAXBuilder;
import org.jdom2.input.StAXEventBuilder;
import org.jdom2.input.sax.XMLReaders;
import org.jdom2.xpath.XPathExpression;
import org.jdom2.xpath.XPathFactory;


public class ReadArrayExpress {

	
	
	public static void main(String[] args) throws IOException, Exception {

		
		
		File xmlFile = new File("/home/mahogny/umeå/project/bias/input/history/experiments.xml");
		
		
		SAXBuilder builder = new SAXBuilder(XMLReaders.NONVALIDATING);
        builder.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd", false);


	    
        Document document = (Document) builder.build(new FileInputStream(xmlFile));

        Element eRoot = document.getRootElement();

		PrintWriter pw=new PrintWriter(new File("/home/mahogny/umeå/project/bias/input/history/ae.csv"));   
		pw.println("id\tdate\texptype");

        for(Element eExperiment:eRoot.getChildren()) {
        	
        	String id=eExperiment.getChildText("id");
        	String releaseDate=eExperiment.getChildText("releasedate");
        	String experimentType=eExperiment.getChildText("experimenttype");
        	
        	pw.println(id+"\t"+releaseDate+"\t"+experimentType);
        }

		pw.close();
		System.out.println("done");
		
	}
}
