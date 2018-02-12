class Score implements Comparable<Score>{
	private double score;
	private boolean valid;
	public Score(double score, boolean valid){
		this.score = score;
		this.valid = valid;
	}
	
	public int compareTo(Score s){
		if (this.score > s.getScore())
			return 1;
		if (this.score == s.getScore())
			return 0;
		return -1;
	}
	
	public double getScore(){return this.score;}
	
	public boolean getValid(){return this.valid;}
	
}
