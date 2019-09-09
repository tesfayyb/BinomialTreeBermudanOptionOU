
public class BinomialTreeBermudanOption {

	
//	implementation of the binomial approximation of the value of a put option on stock for the
//	Ornstein-Uhlenbeck (OU) process.

	/*
	 * S_0 = initial stock price
	 * mu = long term mean stock price
	 * r = risk free interest rate
	 * sigma = Black-Scholes volatility
	 * beta = rate of mean reversion
	 * K = strike price
	 * n = number of time steps
	 * m = number of Bermudan exercise times
	 */
	
	
	
	public static double computeBermudanOption (double S_0, double mu, double sigma, double r,double beta,double K, double T, int n, int m) 
	{
		double dt = T / n ;    //size of time step
		double q ;            //risk neutral probability of an up move (risk neutral measure)
		double s ;           //stock price at time t
		double sigma_OU = sigma * mu ; // Ornstein-Uhlenbeck(OU) volatility
		int i ;              // number of downward stock price movements since time zero
		int j ;              //number of time steps on the binomial tree
		
		double [] bermudanPayoffs = new double [n + 1];  // array to store Bermudan put option payoffs
		
		
		// compute stock prices and put payoffs at time T
		
		for(i = 0; i <= n; i++) {
			
			s = S_0 + n  * beta * (mu - S_0) * dt + (n - 2 * i) * (sigma_OU * Math.sqrt(dt) );
			
			if (s < 0) s = 0;   // stock price can not be negative
			
			bermudanPayoffs[i] = Math.max(K - s, 0.0);
		}
		
		// backward recursion through the tree
		
		for(j = n - 1; j >= 0; j--) {
			for(i = 0; i <= j; i++) {
				
				s = S_0 + j  * beta * (mu - S_0) * dt + (j - 2 * i) * (sigma_OU * Math.sqrt(dt) );
				if (s < 0) s = 0;
				
				q = 0.5 + beta * ((mu - s) -  (mu - S_0)) * Math.sqrt(dt) / (2 * sigma_OU); //risk neutral measure
				
				if (q < 0) q = 0;  // q can not be negative
				if (q > 1) q = 1;  // q can not be greater than 1
				
				if((j)%(n/m) == 0)
				{
					bermudanPayoffs[i] = Math.max(K - s, Math.exp(-r * dt) * (q * bermudanPayoffs[i] + (1 - q) * bermudanPayoffs[i + 1]));
				}
				else
				{
					bermudanPayoffs[i] = Math.exp(-r * dt) * (q * bermudanPayoffs[i] + (1 - q) * bermudanPayoffs[i + 1]);
				}
				
				
			}
		}
		
		
		return bermudanPayoffs[0];
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		double bermudanPutPrice = computeBermudanOption(45, 50, 0.3, 0.05, 0.5, 50, 1, 1000, 4);
		
		System.out.println( "The price of Bermudan put option is :   " + bermudanPutPrice);

	}

}
