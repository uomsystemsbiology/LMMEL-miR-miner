package infodynamics.measures.continuous;

/**
 * <p>This specifies the interface for implementations of
 *  the transfer entropy (see Schreiber, PRL, 2000)
 *  and local transfer entropy (see Lizier et al, PRE, 2008).
 *  This is the <i>apparent</i> transfer entropy; i.e.
 *  we compute the transfer that appears to come from a single
 *  source variable, without examining any other potential sources
 *  (see Lizier et al, PRE, 2008).</p>
 * 
 * <p>Specifically, this specifies the interface for computing
 * the transfer entropy for <i>continuous</i>-valued variables.</p>
 * 
 * @see "Schreiber, Physical Review Letters 85 (2) pp.461-464, 2000;
 *  <a href='http://dx.doi.org/10.1103/PhysRevLett.85.461'>download</a>
 *  (for definition of transfer entropy)"
 * @see "Lizier, Prokopenko and Zomaya, Physical Review E 77, 026110, 2008;
 * <a href='http://dx.doi.org/10.1103/PhysRevE.77.026110'>download</a>
 *  (for definition of <i>local</i> transfer entropy and qualification
 *  of naming it as <i>apparent</i> transfer entropy)"
 *  
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 *
 */
public interface TransferEntropyCalculator extends ChannelCalculator {

	/**
	 * Property name to specify the history length k
	 */
	public static final String K_PROP_NAME = "k_HISTORY";
	
	/**
	 * Initialise the calculator for re-use with new observations.
	 * A new history length k can be supplied here; all other parameters
	 * remain unchanged.
	 * 
	 * @param k history length to be considered.
	 * @throws Exception
	 */
	public void initialise(int k) throws Exception;
	
}
