package genepub;

public class Edge {

	public Edge(String from2, String to2, Double double1) {
		this.from=from2;
		this.to=to2;
		this.pident=double1;
	}
	
	public String from, to;
	public double pident;
	
	@Override
	public boolean equals(Object obj) {
		if(obj instanceof Edge) {
			Edge b=(Edge)obj;
			return (from.equals(b.from) && to.equals(b.to)) ||
					(from.equals(b.to) && to.equals(b.from));
		}
		return super.equals(obj);
	}
	
	@Override
	public int hashCode() {
		return from.hashCode()+to.hashCode();
	}
}
