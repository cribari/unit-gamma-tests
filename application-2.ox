#include<oxstd.oxh>
#import<maximize>

/* Global variables */

static decl s_vY;	// response, dependent variable 
static decl s_tX; 
static decl s_rH;	
static decl s_tXl;	
static decl s_rHl;	
static decl nobs;	 // number of observations
static decl p;	     // number of coefficients in the \mu submodel
static decl q;	     //  number of coefficients in the \phi submodel

const decl dados = "data_application-1.mat";  // data 

const decl l = 1;	// number of restrictions (\phi submodel) 

const decl  beta0 = <0>;   // value of the delta's under H0

/* Log-likelihood function (unrestricted) */

floglik(const vtheta, const adFunc, const avScore, const amHess){
	decl eta = s_tX*vtheta[0:(p-1)];
	decl mu = exp(eta) ./ (1.0 + exp(eta));
	decl zeta =  s_rH*vtheta[p:(p+q-1)] ;
	decl phi = exp(zeta);
	decl ys = log(s_vY);
	decl d = (mu.^(1.0./phi)) ./ (1.0 - (mu.^(1.0./phi)));
	decl ms = -phi ./ d;
	decl yd = log(-log(s_vY));
	decl md = polygamma(phi,0) - log(d);
	decl T1 = diag(exp(eta) ./( (1.0 + exp(eta)) .^2));  // logit link for \mu 
	decl T2 = diag(phi);	        // log link for \phi 
	decl PI = diag(1.0./phi);
	decl M = diag(mu);
	decl MI = diag(1.0 ./ mu);
	decl D = diag(d);
	decl I = unit(nobs,nobs);
	decl P = diag(((-d.*log(mu)) ./ (phi.^2) ) .* (1.0 +d));
  
	adFunc[0] = double(sumc(phi.*log(d) - loggamma(phi) + (d - 1.0) .* ys + (phi - 1.0) .* yd)); 
                             
	if(avScore)
	{
		(avScore[0])[0:(p-1)] = s_tX'T1*PI*MI*D*(I+D)*(ys-ms);
		(avScore[0])[p:(p+q-1)] = s_rH'*T2*(P*(ys-ms)+(yd-md));
	}

	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) )
		return 0;  // failure
	else
		return 1; // sucess
}

/* Log-likelihood function (restricted) */

floglikl(const vtheta, const adFunc, const avScore, const amHess){
    decl p=columns(s_tXl);
    decl q=columns(s_rH);
	decl eta = s_tXl*(vtheta[0:(p-1)]);
	decl mu = exp(eta)./ (1.0 + exp(eta));
	decl zeta =  s_rH*(vtheta[p:(p+q-1)]) ;  
	decl phi = exp(zeta);
	decl ys = log(s_vY);
	decl d = (mu.^(1.0./phi))./(1.0 - (mu.^(1.0./phi)));
	decl ms = -phi./ d;
	decl yd = log(-log(s_vY));
	decl md = polygamma(phi,0) - log(d);
	decl T1 = diag(mu .* (1-mu));
	decl T2 = diag(phi);
	decl PI = diag(1.0./phi);
	decl M = diag(mu);	                                                
	decl MI = diag(1.0./mu);
	decl D = diag(d);
	decl I = unit(nobs,nobs);
	decl P = diag(((-d.*log(mu))./phi.^2 ) .* (1.0+d));
  
	adFunc[0] = double(sumc(phi.*log(d) - loggamma(phi) + (d - 1).*ys + (phi - 1).*yd));
                             
	if(avScore)
	{
		(avScore[0])[0:(p-1)] = s_tXl'T1*PI*MI*D*(I+D)*(ys-ms);
		(avScore[0])[(p):(p+q-1)] = s_rH'*T2*(P*(ys-ms)+(yd-md));		 
	}

	if( isnan(adFunc[0]) || isdotinf(adFunc[0]))
	    return 0; // failure
	else
	    return 1; // success
}

/* Log-lihelihood function that only includes the intercept (used to compute the
   pseudo-R^2 */
   
floglikr(const vtheta, const adFunc, const avScore, const amHess){	
	decl p=columns(s_tXl);
	decl q=columns(s_rHl); 
	decl eta = s_tXl*vtheta[0:(p-1)];
	decl zeta = s_rHl*vtheta[p:((p+q)-1)];
	decl mu = exp(eta) ./ (1.0+exp(eta));
	decl phi = exp(zeta);
	decl ys = log(s_vY);
	decl d = (mu.^(1.0./phi))./(1.0 - (mu.^(1.0./phi)));
	decl ms = -phi./ d;
	decl yd = log(-log(s_vY));
	decl md = polygamma(phi,0) - log(d);
	decl T1 = diag(mu .* (1-mu));
	decl T2 = diag(phi);
	decl PI = diag(1.0./phi);
	decl M = diag(mu);	                                                
	decl MI = diag(1.0./mu);
	decl D = diag(d);
	decl I = unit(nobs,nobs);
	decl P = diag(((-d.*log(mu))./phi.^2 ) .* (1.0+d));

	adFunc[0] = double(sumc(phi.*log(d) - loggamma(phi) + (d - 1).*ys + (phi - 1).*yd));

	if(avScore)
	{
		(avScore[0])[0:(p-1)] = s_tXl'T1*PI*MI*D*(I+D)*(ys-ms);
		(avScore[0])[p:((p+q)-1)] = s_rHl'*T2*(P*(ys-ms)+(yd-md));
	}

	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) )
		return 0;  // failure
	else
		return 1; // success
}
					
nuisance(const matriz){
    decl mat = matriz;
    mat = mat[][:p-l-1]~mat[][p:p+q-1];
    mat = mat[:p-l-1][]|mat[p:p+q-1][];
    return mat;
}

/* main : This is the executable code */ 
 
main(){
	decl data;
	decl dExecTime, controle, ci, i, j, exec=0, contador,  auxiliar;
	decl beta, delta, theta, eta, zeta, mu, phi;
	decl vtheta0, converge, dfunc1, vtheta1,convergel, dfunc2, ys, yd, ys1;
	decl betaols, betaolsl, etaols, muols, varols, phiols,deltaols, etaolsl, muolsl, varolsl, phiolsl;
	decl LR, Skov1,Skov2;
	decl pseudoR2LR, pseudoR2LRc,converger,vp, dfuncl ;
	
	dExecTime = timer();	

	data = loadmat(dados);
	decl dislex = data[][2]; decl QI = data[][3]; 
	s_tX = 1~dislex~QI.^2~dislex.*(QI.^2);  // mean submodel 
	s_rH = 1~dislex~QI; // precision submodel (GammaModel -- II)
  
	s_vY = data[][1];
	nobs = rows(data);
	p = columns(s_tX);
	q = columns(s_rH);

	if(max(s_vY)>= 1.0 || min(s_vY) <= 0.0){
		println("\n\nError: Data outside the standard unit interval (0,1).\n\n");
		exit(2);
	}

	if(p >= nobs){
		println("\n\nError: Number of covariates cannot exceed sample size.\n\n");
		exit(3);
	}

	s_tXl =  s_tX[][0:(p-l-1)];	// matrix X (restricted) 
	
	ys = log(s_vY);
	yd = log(-log(s_vY));
	ys1 = log(s_vY ./ (1.0 - s_vY));

	if(p>1){
		ols2c(ys1, s_tX, &betaols);
	}
	else if(p==1){
		betaols = meanc(ys1);
	}
	
	etaols = s_tX*betaols;
	muols = exp(etaols) ./ (1.0 + exp(etaols));
	varols = ((ys1-etaols)'*(ys1-etaols)) ./
	              ((nobs-p) .* ((1 ./(muols .*(1-muols))) .^2));
	phiols = (muols .* (1.0 - muols) ./ varols) - 1.0;

	if(q>1){
		ols2c(log(phiols), s_rH, &deltaols);
	}
	else if(q==1){
		deltaols = meanc(phiols);
	}
	   	 
	vtheta0 = betaols | deltaols;  // starting values 

	// unrestricted maximum likelihood estimation
	converge = MaxBFGS(floglik, &vtheta0, &dfunc1, 0, FALSE);

	if(converge != MAX_CONV && converge != MAX_WEAK_CONV){
		print("\nConvergence failure!"); 
	}

	decl etahat = s_tX * vtheta0[0:(p-1)];
	decl muhat = exp(etahat) ./ (1+exp(etahat));
	decl T1 = diag(muhat .* (1-muhat));	 // 1/g'(mu)
	decl zetahat = s_rH * vtheta0[p:(p+q-1)];
	decl phihat = exp(zetahat);
	decl T2 = diag(phihat);
	decl s_tys = diag(ys);
	decl s_tyd = diag(yd);
	decl M = diag(muhat);
	decl IM = diag(1.0./muhat);
	decl F = diag(phihat);
	decl IF = diag(1.0 ./ phihat);
	decl dhat = (muhat.^(1.0./phihat))./(1.0 - ( muhat.^(1.0./phihat)));
	decl Dhat = diag(dhat);
	decl MS = diag(-phihat ./ dhat);
	decl VS = diag(phihat ./ (dhat.^2));
	decl MD = diag(polygamma(phihat,0) - log(dhat));
	decl VD = diag(polygamma(phihat,1));
	decl C = diag(-1.0 ./ dhat);
	decl S1 = diag((2.0*muhat - 1.0)./((muhat .*(1-muhat)) .^2));
	decl S2 = diag(-1.0 ./ (phihat .^2));
	decl Phat = diag(((-dhat.*log(muhat))./(phihat.^2) .* (1+dhat)));
	decl K1 = diag(dhat .* (1+dhat)./(phihat .* muhat.^2) .* ((1.0+2.0.*dhat)./phihat - 1.0));
    decl K2 = diag(dhat .* (1+dhat)./(phihat .* muhat.^2 .* muhat.^(1.0 ./ phihat)));
	decl K3 = diag(dhat.*(1+dhat) ./ (muhat.* phihat));
	decl Z1 = diag((-dhat.*(1+dhat) ./( muhat .* phihat) .* log(muhat) .* ((1.0 + 2.0 .*dhat)./(phihat.^2))) - (dhat.*(1.0+dhat)./(muhat.*phihat.^2)));
	decl Z2 = diag(dhat.*(1+dhat) ./( muhat .* phihat) .* (1.0 ./ dhat + (log(muhat)./phihat) + (log(muhat)./(dhat.*phihat)) ));
	decl L1 = diag(dhat ./ (phihat .^3) .* (log(muhat)./phihat + dhat.*log(muhat)./phihat + 2.0) .* log(muhat) .* (1.0+dhat) + ((dhat.^2 .* log(muhat).^2 ./ (phihat.^4)) .* (1.0+dhat))) ;
    decl L2 = diag(-dhat.*(1+dhat).* log(muhat) ./( phihat.^2) .* (1.0 ./ dhat + (log(muhat)./phihat) + (log(muhat)./(dhat.*phihat) )) - polygamma(phihat,1) - log(muhat).*(1.0+dhat)./(phihat.^2));
	decl L3 = diag(-dhat.*log(muhat)./(phihat.^2) .* (1.0+dhat));
	decl Ibetabeta = s_tX'* K2 * T1^2 * s_tX;
	decl Ibetadelta = -s_tX'* Z2 * T1 * T2 * s_rH;
	decl Ideltabeta = Ibetadelta';
	decl Ideltadelta = -s_rH'*T2*L2*T2*s_rH;
	decl fisher = (Ibetabeta~Ibetadelta)|(Ideltabeta~Ideltadelta);
	decl Ifisher = invert(fisher);	  
	decl Jbetabeta = -s_tX'*((K1*(s_tys - MS) - K2)*T1 - S1*T1^2*K3*(s_tys - MS))*T1*s_tX;
	decl Jbetadelta = -s_tX'*(Z1*(s_tys - MS) + Z2)*T1*T2*s_rH;
	decl Jdeltabeta = Jbetadelta';
	decl Jdeltadelta = -s_rH'*(T2*(L1*(s_tys - MS)+L2) - (L3*(s_tys - MS)+(s_tyd-MD)) *S2*T2^2)*T2*s_rH;
	decl infobs = (Jbetabeta~Jbetadelta)|(Jdeltabeta~Jdeltadelta);
	decl infobsinv = invert(infobs);
	decl stderrors = sqrt(diagonal(Ifisher));	 

	p = columns(s_tXl);
	if(p>1){
		ols2c(ys1, s_tXl, &betaolsl);
	}
	else if(p==1){
		betaolsl = meanc(ys1);
	}
	etaolsl = s_tXl*betaolsl;
	muolsl = exp(etaolsl) ./ (1.0 + exp(etaolsl));
	varolsl = ((ys1-etaolsl)'*(ys1-etaolsl)) ./
	              ((nobs-p) .* ((1 ./(muolsl .*(1-muolsl))) .^2));
	phiolsl = (muolsl .* (1.0 - muolsl) ./ varolsl) - 1.0;

	if(q>1){
		ols2c(log(phiols), s_rH, &deltaols);
	}
	else if(q==1){
		deltaols = meanc(phiols);
	}
	   	 
	vtheta1 =  betaolsl | deltaols;  // starting values 
	
	// restricted maximum likelihood estimation
	convergel = MaxBFGS(floglikl, &vtheta1, &dfunc2, 0, FALSE);
	   
	if(convergel != MAX_CONV && convergel != MAX_WEAK_CONV){
		print("\nConvergence failure!");
	}
	
	decl etatil = s_tX*(vtheta1[0:(p-1)]|beta0');	
	decl mutil = exp(etatil) ./ (1+exp(etatil));
	decl T1til = diag(mutil .* (1-mutil));
	decl zetatil = s_rH*(vtheta1[p:(p+q-1)]);
	decl phitil = exp(zetatil);
	decl T2til = diag(phitil);
	decl Mtil = diag(mutil);
	decl IMtil = diag(1.0 ./ mutil);
	decl Ftil = diag(phitil);
	decl IFtil = diag(1.0 ./ phitil);
	decl dtil = (mutil.^(1.0 ./ phitil))  ./ (1.0 - (mutil.^(1.0 ./ phitil)));
	decl Dtil = diag(dtil);
	decl mstil = -phitil ./ dtil;
	decl MStil = diag(mstil);
	decl mdtil = polygamma(phitil,0) - log(dtil);
	decl MDtil = diag(mdtil);
	decl VStil = diag(phitil ./ (dtil.^2));
	decl VDtil = diag(polygamma(phitil,1));
	decl Ctil = diag(-1.0 ./ dtil);
	decl S1til = diag((2*mutil-1) ./ ((mutil .* (1-mutil)).^2));
	decl S2til = diag(-1.0 ./ (phitil .^2));
	decl K1til = diag(dtil .* (1+dtil)./(phitil .* mutil.^2) .* ((1+2.0.*dtil)./phitil - 1.0));
	decl K2til = diag(dtil .* (1+dtil)./(phitil .* mutil.^2 .* mutil.^(1.0 ./ phitil)));
	decl K3til = diag(dtil.*(1+dtil) ./ (mutil .* phitil));
	decl Z1til = diag(-dtil.*(1+dtil) ./( mutil .* phitil) .* log(mutil) .* ((1.0 + 2.0 .*dtil)./(phitil.^2)) - dtil.*(1+dtil)./(mutil .* phitil.^2));
	decl Z2til = diag(dtil.*(1+dtil) ./( mutil .* phitil) .* (1.0 ./ dtil + (log(mutil)./phitil) + (log(mutil)./(dtil.*phitil)) ));
	decl L1til = diag(dtil ./ (phitil .^3) .* (log(mutil)./phitil + ( dtil.*log(mutil)./phitil) + 2) .* log(mutil) .* (1+dtil) + (dtil.^2 .* log(mutil).^2 ./ (phitil.^4)) .* (1+dtil)) ;
	decl L2til = diag(-dtil.*(1+dtil).* log(mutil) ./( phitil.^2) .* (1.0 ./ dtil + (log(mutil)./phitil) + (log(mutil)./(dtil.*phitil) )) - polygamma(phitil,1) - log(mutil).*(1+dtil)./phitil.^2);
	decl L3til = diag(-dtil.*log(mutil)./(phitil.^2) .* (1+dtil));
	decl Ibetabetatil = s_tX'* K2til * T1til^2 * s_tX;
	decl Ibetadeltatil = -s_tX'* Z2til * T1til * T2til * s_rH;
	decl Ideltabetatil = Ibetadelta';
	decl Ideltadeltatil = -s_rH'*T2til*L2til*T2til*s_rH;
	decl fishertil = (Ibetabetatil~Ibetadeltatil)|(Ideltabetatil~Ideltadeltatil);
	decl Ifishertil = invert(fishertil);
	decl Jbetabetatil = -s_tX'*((K1til*(s_tys - MStil) - K2til)*T1til - S1til*T1til^2*K3til*(s_tys - MStil))*T1til*s_tX;
	decl Jbetadeltatil = -s_tX'*(Z1til .*(s_tys - MStil) + Z2til)*T1til*T2til*s_rH;
	decl Jdeltabetatil = Jbetadelta';
	decl Jdeltadeltatil = -s_rH'*(T2til*(L1til*(s_tys - MStil)+L2til) - (L3til*(s_tys - MStil)+(s_tyd-MDtil)) *S2til*T2til^2)*T2til*s_rH;
	decl infobstil = (Jbetabetatil~Jbetadeltatil)|(Jdeltabetatil~Jdeltadeltatil);
	decl infobsinvtil = invert(infobstil);
	decl Iden = unit(nobs,nobs);
	decl Ptil = diag(((-dtil.*log(mutil))./phitil.^2 ) .* (1.0 + dtil));
	decl escorebetatil = s_tX'*T1til*IFtil*IMtil*Dtil*(Iden+Dtil)*(ys-mstil);
	decl escoredeltatil = s_rH'*T2til*(Ptil*(ys-mstil) + (yd-mdtil));
	decl escoretil = escorebetatil | escoredeltatil;

	//****************** pseudo-R^2 (likelihod-based) 

	s_tXl = ones(nobs, 1);
	s_rHl = ones(nobs, 1);
	decl ynew = log( s_vY ./ (1.0-s_vY) );
	decl ynewbal = meanc(ynew);
	decl muhatl = exp(ynewbal)/(1.0 + exp(ynewbal));
	decl phi_l = ((1.0/(varc(ynew)*muhatl*(1.0-muhatl)))); 
	decl v_deltal =(log(phi_l));
	 
	vp = ynewbal|v_deltal;        
	p = columns(s_tX);
	
	/* restricted maximum likelihood estimation */
 	converger = MaxBFGS(floglikr, &vp, &dfuncl, 0, FALSE);                                                           
 	                                                                                                         
 	if(converger == MAX_CONV || MAX_WEAK_CONV){                                                                                                       
 		pseudoR2LR = 1.0 - (exp(dfuncl)/exp(dfunc1))^(2/nobs);                                                
		pseudoR2LRc = 1.0 - (1 - pseudoR2LR)*((nobs - 1)/(nobs - (1 + 0.4)*(p + 1) - (1 - 0.4)*(q + 1)));    
	}
	
    //******************** pseudo-R^2 (Ferrari and Cribari-Neto)
   
	decl g_liber = nobs-(p+q);
	decl pseudoR2 = (correlation(ynew~etahat)[0][1])^2;
	decl pseudoR2c = 1.0 - (1 - pseudoR2)*((nobs - 1)/g_liber);  			   
			   
	//******************** Model selection criteria
	 
    decl AIC = double(-2.0*dfunc1 + 2.0*(p+q))	;
	decl BIC = double(-2.0*dfunc1 + ((p+q)*log(nobs)));
	decl AICc = double(-2.0*dfunc1 + 2.0*(p+q)*(nobs/(nobs-(p+q)-1)));

	//********************* Quantities used to compute the modified test statistics
	
	decl iota = ones(nobs,1);
	decl q1 = s_tX'*T1*IF*IM*Dhat*(Iden+Dhat)*(VS*(Dhat-Dtil) + C*(F-Ftil))*iota;
	decl q2 = s_rH'*T2*((Phat*VS + C)*(Dhat-Dtil) + (Phat*C + VD)*(F-Ftil))*iota;
	decl qb = q1|q2;
	decl upsilon1 = s_tX'*T1*IF*IM*Dhat*(Iden+Dhat)*VS*(Iden+Dtil)*Dtil*IMtil*IFtil*T1til*s_tX;
	decl upsilon2 = s_tX'*T1*IF*IM*Dhat*(Iden+Dhat)*(VS*Ptil+ C)*T2til*s_rH;
	decl upsilon3 = s_rH'*T2*(Phat*VS + C)*(Iden+Dtil)*Dtil*IMtil*IFtil*T1til*s_tX;
	decl upsilon4 = s_rH'*T2*(Phat*VS*Ptil +(Phat+Ptil)*C + VD)*T2til*s_rH;
	decl upsilon = (upsilon1~upsilon2)|(upsilon3~upsilon4);
	decl upsiloninv = invert(upsilon);
	decl nuiinfobstil = nuisance(infobstil);
	decl nuisance2 = nuisance(fishertil*upsiloninv*infobs*Ifisher*upsilon);

	// Test statistics: LR, Skov1, Skov2
	LR = 2*(dfunc1-dfunc2);	// LR

	decl xi = fabs(fabs(determinant(fishertil))^(1/2)*fabs(determinant(fisher))^(1/2)*
	      		fabs(determinant(upsilon))^(-1)*fabs(determinant(nuiinfobstil))^(1/2)*
				fabs(determinant(nuisance2))^(-1/2)*
				fabs(escoretil'*upsiloninv*fisher*infobsinv*upsilon*Ifishertil*escoretil)^(l/2)/
				fabs((LR)^((l/2)-1)*escoretil'*upsiloninv*qb));
								
	Skov1 = LR-2*log(xi);
	Skov2 = LR*(1-(log(xi)/LR))^2;

	println("\nOX PROGRAM: ", oxfilename(0));
	println("\nOX VERSION: ", oxversion() );
	println("\nDISTRIBUTION: Unit Gamma" );
	println("\nFILE THAT CONTAINS THE DATA: ", dados);							
	println("\nSAMPLE SIZE: ", nobs);
	println("\nTOTAL NUMBER OF COVARIATES: ", p+q );
	println("\nConvergence status, unrestricted: ", MaxConvergenceMsg(converge));
	println("\nConvergence status, restricted: ", MaxConvergenceMsg(convergel));

	println("\nNumber of parameters in the submodel for mu (p): ", p);
	println("\nNumber of parameters in the submodel for phi (q): ", q);
	println("\nNumber of restrictions (l): ", l);		  
	
	println("\nNull hypothesis H0: (beta_",(p),") = 0");

	println("\nParameter estimates (beta, delta):" );
	println("\n", "%8.3f", "%r",{"estimate","standard error"}, "%c",{"beta_1","beta_2","beta_3","beta_4","beta_5",
	        "beta_6", "beta_7","beta_8","beta_9","beta_10"},	vtheta0[:(p-1)]'|stderrors[:(p-1)]);
	println("%8.3f",	"%r",{"estimate","standard error"}, "%c",{"delta_1","delta_2","delta_3","delta_4","delta_5",
	        "delta_6", "delta_7","delta_8","delta_9","delta_10"}, vtheta0[p:(p+q-1)]'|stderrors[p:(p+q-1)]);
			  
	decl resul = zeros(2,3);
	resul[0][0] = LR;
	resul[0][1] = Skov1;
	resul[0][2] = Skov2;
	resul[1][0] = 1-probchi(LR, l);
	resul[1][1] = 1-probchi(Skov1, l);
	resul[1][2] = 1-probchi(Skov2, l);

	println("Test statistics: ", "%10.4f", "%c", {"w", "w*", "w**"}, "%r", {"statistic", "p-value"}, resul);

	println("Pseudo-R^2: ", "%10.4f", "%c", {"R2", "R2c", "R2LR", "R2LRc"},
	pseudoR2~pseudoR2c~pseudoR2LR~pseudoR2LRc);

	println("Model selection criteria: ", "%10.4f", "%c", {"AIC", "BIC", "AICc"}, AIC~BIC~AICc);

	println("DATE: ", date());
	println("\nTIME: ", time(),"\n");
	println("\nEXECUTION TIME: ", timespan(dExecTime));
	println("\n");
}