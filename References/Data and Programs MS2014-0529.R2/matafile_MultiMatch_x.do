
/*
august 2: 	this version eliminates  regressors that generate collinearity. 
august 3: 	seed is fixed for bootstrap. 
			nboot is now an option
		

*/
mata:

mata clear
mata set matastrict on


void estimate_multi_match(string scalar varname, string scalar  rfvar, string scalar ivvar, string scalar idvar, string scalar cluvar,
string scalar touse, real scalar pd, real scalar db, 
real  scalar K1, real  scalar K2, real scalar K3, real scalar modeltype, real scalar ak, real scalar wls,
real scalar bs, real scalar cludum)
/* this will eventually become the estimation routine...  */

{
	real matrix  x1, x2, x3, bodat, bodattemp, nn, vhat,idel, cludummies
	real vector y, id, bhat, KK,ideleted, clux, cluvals
	real scalar n, fopt, n1, nx, nper, nclu, i, j, K4
	string vector names
	string scalar allvars

/*	rseed(1771546426543) */
	

	allvars= varname +" "+  rfvar+" " +ivvar+ " " + idvar + " " + cluvar 
	

	
	
	/* allvars */

	st_view(bodattemp=.,.,tokens(allvars),touse)
	
	bodat=bodattemp
	
	K4=0
	if (cludum==1) {
		cluvals=uniqrows(bodat[.,cols(bodat)]) 
		nclu=rows(cluvals)
		cludummies=J(rows(bodat),nclu-1,0)
		for (i=1;i<nclu; i++){
			for (j=1;j<=rows(bodat); j++){
				if(bodat[j,cols(bodat)]==cluvals[i]) {
					cludummies[j,i]=1
				}
				
			}
		}
		K4=nclu-1
	}
	
	if (K4==0) {
		bodat=bodat[.,1],bodat[.,K2+2..K2+K1+1],bodat[.,2..K2+1],bodat[.,(K2+K1+2)..cols(bodat)]
	}
	else {
		bodat=bodat[.,1],bodat[.,K2+2..K2+K1+1],bodat[.,2..K2+1],cludummies,bodat[.,(K2+K1+2)..(cols(bodat)-2)],cludummies, bodat[.,(cols(bodat)-1)..cols(bodat)]
		K2=K2+K4
		K3=K3+K4
	}
		
	/* Note that the last column is cluster. the one before that is id */
	
	
	if (pd==1){
	/* pd==1 means that we want to print out details for debugging */ 
		"MODELTYPE:"
		modeltype
	}
	if (modeltype==0){
	/* some of the first part is redundant (but hopefully harmless) it was used to catch mistakes in earlier versions */
	/* The following is wrong if K1=0. But that is ruled out in the calling routine */
		
		if (K1==0) {
			K1=(cols(bodat)-3-K2)
		}
		else if (K2==0) {
			K2=cols(bodat)-3-K1
		}
		
		K3=K2
		
		bodat=bodat[.,1..1+K1], bodat[.,2+K1..1+K1+K2], bodat[.,2+K1..1+K1+K2], bodat[.,(cols(bodat)-1)..cols(bodat)]
		/* the next line just makes regression a special case of IV */
		modeltype=2
	}
	
	if (modeltype==1) {
	/* some of the first part is redundant (but hopefully harmless) it was used to catch mistakes in earlier versions */
	/* The following is wrong if K1=0. But that is ruled out in the calling routine */
		K3=0
		if (K1==0) {
			K1=cols(bodat)-3-K2
		}
		else if (K2==0) {
			K2=cols(bodat)-3-K1
		}
		nx=K2
		
		 if  ((K1+K2+K3+3)!=cols(bodat)){
			"K1,K2 or K3 is not right 001" 
		 }
		 else{
			 KK=(K1,K2,K3,K4)
			 cleanx2(bodat, KK, ideleted,pd)
			 estimate_logit(bodat,KK, pd, db, ak, wls, bs, nper, bhat, vhat)
			 
		 }
	}
	else if (modeltype==2) {
		if (K1==0) {
			K1=cols(bodat)-3-K2-K3
		}
		else if (K2==0) {
			K2=cols(bodat)-3-K1-K3
		}
		else if (K3==0) {
			K3=cols(bodat)-3-K1-K2
		}
		nx=K2
		 if  ((K1+K2+K3+3)!=cols(bodat)){
			"K1,K2 or K3 is not right 002" 
		 }
		 else{
		 KK=(K1,K2,K3,K4)
		 if (K3<K2) {
		 "ERROR: Underidentified"
		 }
		 else{
			ideleted=0
			cleanx2(bodat, KK, ideleted,pd)	
			estimate_2sls(bodat,KK, pd, db, ak, wls, bs, nper, bhat, vhat)
			 
		 }
		 }
	}
	else {
	"ERROR: modeltype not implemented"
	}
	
	idel=ideleted'
	 
	bhat=bhat[1..(cols(bhat)-K4)]
	vhat=vhat[1..(rows(vhat)-K4),1..(cols(vhat)-K4)]
	 
	
	st_matrix("r(beta)", bhat)
	st_matrix("r(V)", vhat)
	st_matrix("r(N)", n\nper\modeltype\fopt\colsum(ideleted'))
	st_matrix("r(K)", (KK,nx))
	st_matrix("r(I)", idel)
	
}

void bo_organizedata(real matrix bodat, real vector KK, real scalar db,
		real matrix nn, real vector y, real matrix x1, real matrix x2, 
		real matrix x3, real vector id, real vector clu)
{
	
	real scalar i,j,n, K1, K2, K3
	
	K1=KK[1]
	K2=KK[2]
	K3=KK[3]

	n=rows(bodat)
	if (cols(bodat)!=(1+K1+K2+K3+1+1)){
		"colums dont match"
	}
	else{
		if (db==1){
			"....so far so good... colums match"
		}
	}
	if (cols(nn)==1){
		bodat=bodat,(1..rows(bodat))' /* I dont quite remember what this does. But why not ;-) */
		bodat=sort(bodat,(cols(bodat)-2,cols(bodat)-1,cols(bodat))) 
		bodat=bodat[.,1..cols(bodat)-1]
		i=2
		j=1
		nn=J(rows(bodat),2,0)
		nn[1,1]=1
		for (i=2; i<=rows(bodat); i++) {
			if (bodat[i,cols(bodat)-1] != bodat[i-1,cols(bodat)-1]) {
				nn[j,2]=i-1
				j=j+1
				nn[j,1]=i
			}
		}
		nn[j,2]=rows(bodat)
		nn=nn[1..j,.]
		nn=(nn,(nn[.,2]-nn[.,1]:+1))
	}
	/* 
	The data has now been sorted and nn is a table of starting and ending rows numbers for each invididual. 
	The third colum is the number of observations of that individual.
	*/	
	y=bodat[.,1]
	n=rows(y)
	x1=J(n,1,1)
	if (K1>0){
		x1=x1,bodat[.,2..(K1+1)]
	}
	x2=J(n,1,1)
	if (K2>0){
		x2=x2,bodat[.,(K1+2)..(K1+1+K2)]
	}
	x3=J(n,1,1)
	if (K3>0){
		x3=x3,bodat[.,(K1+K2+2)..(K1+1+K2+K3)]
	}
	id=bodat[.,cols(bodat)-1]
	clu=bodat[.,cols(bodat)]

}

void cleanx2(real matrix bodat, real vector K, real vector ideleted, real scalar pd)
/* takes bodat and eliminates colums of x2 that generate colinearity 
ideleted is a dummy fo whether a colums of bodat was deleted */
{	
	real matrix x1,x2, x3, xx2, tmp
	real vector y, id, sv, clu
	real scalar i
	
	y=bodat[.,1]
	x1=bodat[.,2..(K[1]+1)]
	x2=bodat[.,(K[1]+2)..(K[1]+1+K[2])]
	if (K[3]>0){
		x3=bodat[.,(K[1]+K[2]+2)..(K[1]+K[2]+K[3]+1)]
	}
	id=bodat[.,cols(bodat)-1]
	clu=bodat[.,cols(bodat)]
	
	ideleted=J(1,cols(x2),0)
	xx2=J(rows(x2),1,1)
	for (i=1; i<=cols(x2); i++) {
		tmp=xx2,x2[.,i]
		sv=svdsv(tmp) 
		if (sv[1]<1.0d6*sv[cols(tmp)]){
			xx2=tmp
		}
		else {
			if (pd==1){
				"variable i dropped"
				"singular values are"
				sv
			}
			ideleted[i]=1
			K[2]=K[2]-1
		}
	}
	
	bodat=y,x1,xx2[.,2..cols(xx2)]
	if (K[3]>0){
		bodat=bodat,x3
	}
	bodat=bodat,id,clu
	
	
}

void estimate_2sls(real matrix bodat, real vector KK,  
	real scalar pd, real scalar db, real scalar ak, real scalar wls, real scalar nboot,
	real scalar nper, real vector bhat, real matrix vhat)
{
	real matrix bboot , nn
	
	nn=1
	est2sls(bodat,KK, nn, pd, db, ak, wls, bhat)
	if (pd==1){
		"point estimates"
		bhat
	}
	nper=rows(nn)
	boots_2sls( bodat, KK, nn, nboot, bboot, wls, bhat, pd, db, ak)
	vhat=variance(bboot) 
	
}

void est2sls(real matrix bodat, real vector KK, real matrix nn,
	real scalar pd, real scalar db, real scalar ak, real scalar wls,
	real vector bhat)
{
	real matrix x1,x2,x3, gh, s2, xx,zz
	real vector y,id, bh1, bh2, yhat, xhat,yy, e2, Li, sihat, iclu
	real scalar i, epsilo	

	bo_organizedata(bodat, KK, db, nn, y, x1, x2, x3, id, iclu)
	bo_OLS(y, x1, bh1, s2, yhat)	
	bo_OLS(x2, x3, gh, s2, xhat)
	
	xx=x2[nn[.,1],.]
	zz=xhat[nn[.,1],.]
	yy=J(rows(nn),1,0)
	for (i=1; i<=rows(nn); i++){
		yy[i]=colsum(y[nn[i,1]..nn[i,2]])
	}
	
	yy=yy-(nn[.,3]:-1):*yhat[nn[.,1],.]
	 
	bh2=luinv(zz'*xx)*(zz'*yy)
	
	if (wls==1){
		e2=(yy-xx*bh2):^2
		Li=(J(rows(nn),1,1),nn[.,3])
		bo_OLS(e2, Li, bh2, s2, sihat)
		epsilo=1.0d-5
		sihat=1:/sqrt(abs(sihat):+(epsilo*variance(yy)))
		zz=(sihat*J(1,cols(zz),1)):*zz
		bh2=luinv(zz'*xx)*(zz'*yy)
	}
	
	bhat= bh1,bh2'
 
}

void boots_2sls(real matrix bodatoriginal,  real vector KK,  real matrix nnd, 
                   real scalar nboot, real matrix bboot, real scalar wls, real vector bhat, real scalar pd, real scalar db, real scalar ak)
{

	real matrix  tmpdat, bodat, nn, x1,x2,x3, nnclu, bodatd, cludummies
	real vector  iii, y, id, p,idold,cluvals
	real scalar  k, nper, ib, i, j1, j2,j3, j,nclu
	
	if (db==1){
		"WE ARE NOW IN THE BOOTSTRAP ROUTINE"
	}

	bodatd=bodatoriginal
	if (db==1){
		"original"
		bodatd[1..35,(cols(bodatd)-4)..cols(bodatd)]
	}
	bodatd=sort(bodatd,(cols(bodatd)-1)) 
	if (db==1){
		"sorted"
		bodatd[1..35,(cols(bodatd)-4)..cols(bodatd)]
	}
	idold=bodatd[.,(cols(bodatd)-1)]
	j=1
	bodatd[1,(cols(bodatd)-1)]=j
	for (i=2; i<=rows(bodatd); i++) {
		if (idold[i] != idold[i-1]) {
			j=j+1
		}
		bodatd[i,(cols(bodatd)-1)]=j
	}
	/* we have now relabeled id to run from 1,2,3,...*/
	if (db==1){
		"new id"
		bodatd[1..35,(cols(bodatd)-4)..cols(bodatd)]
	}

	bodatd=bodatd,(1..rows(bodatd))' /* I dont quite remember what this does. But why not ;-) */
	bodatd=sort(bodatd,(cols(bodatd)-1,cols(bodatd)-2,cols(bodatd))) 
	bodatd=bodatd[.,1..cols(bodatd)-1]
	i=2
	j=1
	nnclu=J(rows(bodatd),2,0)
	nnclu[1,1]=1
	for (i=2; i<=rows(bodatd); i++) {
		if (bodatd[i,cols(bodatd)] != bodatd[i-1,cols(bodatd)]) {
			nnclu[j,2]=i-1
			j=j+1
			nnclu[j,1]=i
		}
	}
	nnclu[j,2]=rows(bodatd)
	nnclu=nnclu[1..j,.]
	nnclu=(nnclu,(nnclu[.,2]-nnclu[.,1]:+1))
	
	
 	
	k=cols(bhat)	
	bboot=J(nboot,k,0)
	nper=rows(nnclu) /* nper is number of clusters */

	for (ib=1; ib<=nboot; ib++) {		
			
		if (pd==1) {
			" "
			"BOOTSTRAP NUMBER&"
			ib		
			displayflush()
		}

		bodat=J(nper*colmax(nnclu[.,3]),cols(bodatd),0)
		iii=ceil(uniform(nper,1)*nper)
		j1=1
		
		for (i=1; i<=nper; i++) {
			j2=j1+nnclu[iii[i],3]-1
			bodat[j1..j2,.]=bodatd[nnclu[iii[i],1]..nnclu[iii[i],2],.]
			j3=j2-j1+1
			bodat[j1..j2,cols(bodat)]=J(j3,1,i)
			bodat[j1..j2,cols(bodat)-1]=J(j3,1,i):/(nper+1)+bodat[j1..j2,cols(bodat)-1];  
			j1=j2+1
		}
		nn=0
		bodat=bodat[1..j1-1,.]
		
		if (KK[4]>0){
		
			/* this could be done more simply. But this makes it easier to compare to cludummies above */
			cluvals=uniqrows(bodat[.,cols(bodat)]) 
			nclu=rows(cluvals)
			cludummies=J(rows(bodat),nclu-1,0)
			for (i=1;i<nclu; i++){
				for (j=1;j<=rows(bodat); j++){
					if(bodat[j,cols(bodat)]==cluvals[i]) {
						cludummies[j,i]=1
					}
				}
			}
		 
			bodat[.,(KK[1]+1+KK[2]-nclu+2)..(KK[1]+1+KK[2])]=cludummies
			bodat[.,(KK[1]+1+KK[2]+KK[3]-nclu+2)..(KK[1]+1+KK[2]+KK[3])]=cludummies

		}


		est2sls(bodat,KK, nn, pd, db, ak, wls, p)
		bboot[ib,.]=p
	}	
		
}

void estimate_logit(real matrix bodat, real vector KK,  
	real scalar pd, real scalar db, real scalar ak, real scalar wls, real scalar nboot,
	real scalar nper, real vector bhat, real matrix vhat)
{
	real matrix bboot , nn
	
	estimate_logit_a(bodat,KK, nn, pd, db, ak,  nper, bhat, vhat)
	if (nboot>0) {
		boots_logit( bodat, KK, nn, nboot, bboot, bhat, pd, db, ak)
		vhat=variance(bboot) 
	}
	nper=rows(nn)
	
	
}

void estimate_logit_a(real matrix bodat, real vector KK,  real matrix nn, 
	real scalar pd, real scalar db, real scalar ak, 
	real scalar nper, real vector bhat, real matrix vhat)
{
	real matrix bboot, x1, x2, x3
	real vector y , id, iclu
	
			 nn=1
			 bo_organizedata(bodat, KK, db, nn, y, x1, x2, x3, id,iclu)
			 est_logit(y, x1, x2,  nn, pd, db, ak, bhat, vhat)
			 
	
}


void boots_logit(real matrix bodatoriginal,  real vector KK,  real matrix nnd, 
                   real scalar nboot, real matrix bboot, real vector bhat, real scalar pd, real scalar db, real scalar ak)
{

	real matrix  tmpdat, bodat, nn, x1,x2,x3, nnclu, bodatd, vp, cludummies
	real vector  iii, y, id, p,idold, cluvals
	real scalar  k, nper, ib, i, j1, j2,j3, j, nclu
	
	if (db==1){
		"WE ARE NOW IN THE BOOTSTRAP ROUTINE"
	}

	bodatd=bodatoriginal
	if (db==1){
		"original"
		bodatd[1..35,(cols(bodatd)-4)..cols(bodatd)]
	}
	bodatd=sort(bodatd,(cols(bodatd)-1)) 
	if (db==1){
		"sorted"
		bodatd[1..35,(cols(bodatd)-4)..cols(bodatd)]
	}
	idold=bodatd[.,(cols(bodatd)-1)]
	j=1
	bodatd[1,(cols(bodatd)-1)]=j
	for (i=2; i<=rows(bodatd); i++) {
		if (idold[i] != idold[i-1]) {
			j=j+1
		}
		bodatd[i,(cols(bodatd)-1)]=j
	}
	/* we have now relabeled id to run from 1,2,3,...*/
	if (db==1){
		"new id"
		bodatd[1..35,(cols(bodatd)-4)..cols(bodatd)]
	}

	bodatd=bodatd,(1..rows(bodatd))'  
	bodatd=sort(bodatd,(cols(bodatd)-1,cols(bodatd)-2,cols(bodatd))) 
	bodatd=bodatd[.,1..cols(bodatd)-1]
	i=2
	j=1
	nnclu=J(rows(bodatd),2,0)
	nnclu[1,1]=1
	for (i=2; i<=rows(bodatd); i++) {
		if (bodatd[i,cols(bodatd)] != bodatd[i-1,cols(bodatd)]) {
			nnclu[j,2]=i-1
			j=j+1
			nnclu[j,1]=i
		}
	}
	nnclu[j,2]=rows(bodatd)
	nnclu=nnclu[1..j,.]
	nnclu=(nnclu,(nnclu[.,2]-nnclu[.,1]:+1))
	
	
 	
	k=cols(bhat)	
	bboot=J(nboot,k,0)
	nper=rows(nnclu) /* nper is number of clusters */

	for (ib=1; ib<=nboot; ib++) {		
			
		if (pd==1) {
			" "
			"BOOTSTRAP NUMBER&"
			ib		
			displayflush()
		}

		bodat=J(nper*colmax(nnclu[.,3]),cols(bodatd),0)
		iii=ceil(uniform(nper,1)*nper)
		j1=1
		
		for (i=1; i<=nper; i++) {
			j2=j1+nnclu[iii[i],3]-1
			bodat[j1..j2,.]=bodatd[nnclu[iii[i],1]..nnclu[iii[i],2],.]
			j3=j2-j1+1
			bodat[j1..j2,cols(bodat)]=J(j3,1,i)
			bodat[j1..j2,cols(bodat)-1]=J(j3,1,i):/(nper+1)+bodat[j1..j2,cols(bodat)-1];  
			j1=j2+1
		}
		if (db==1){
		"Bootstrap sample constructed"
		 
		}
		
		nn=0
		bodat=bodat[1..j1-1,.]
		
		if (KK[4]>0){
		
			/* this could be done more simply. But this makes it easier to compare to cludummies above */
			cluvals=uniqrows(bodat[.,cols(bodat)]) 
			nclu=rows(cluvals)
			cludummies=J(rows(bodat),nclu-1,0)
			for (i=1;i<nclu; i++){
				for (j=1;j<=rows(bodat); j++){
					if(bodat[j,cols(bodat)]==cluvals[i]) {
						cludummies[j,i]=1
					}
				}
			}
		 
			bodat[.,(KK[1]+1+KK[2]-nclu+2)..(KK[1]+1+KK[2])]=cludummies
			bodat[.,(KK[1]+1+KK[2]+KK[3]-nclu+2)..(KK[1]+1+KK[2]+KK[3])]=cludummies

		}
		
		estimate_logit_a(bodat,KK, nn, pd, db, ak,  nper, p, vp)
		bboot[ib,.]=p
	}	
		
}


void est_logit(real vector y, real matrix x1, real matrix x2, real matrix nn,
	real scalar pd, real scalar db, real scalar ak, 
	real vector bhat, real matrix vhat)

{
	real matrix  H1,H2, sco1,sco2, V, H,sco, H21
	real vector  phat, bhat1,bhat2,bstart,g, delt, df, bb1, df1, sv
	real scalar  i,todo,f,fopt,  delta, kapstart, f1, conv
	transmorphic S 

		/*if y is not discrete (i.e. 0 or 1) then send an error message */
	for (i=1; i<=rows(y); i++) {
		if (y[i]!=1 & y[i]!=0) {
			"**** ---error: use discrete y (0 or 1)--- ****"
			/*deliberately kill the program */
			bo_kill()
			}
		} 
		
	   bstart=J(1,cols(x1),0)
		S=optimize_init()
		if (pd==0) {
			optimize_init_verbose(S,0)
			optimize_init_tracelevel(S, "none")
		}
		optimize_init_conv_maxiter(S,100)
		optimize_init_evaluator(S, &bo_logit())
		optimize_init_evaluatortype(S, "v1")
		optimize_init_technique(S,"nr")
		optimize_init_params(S,bstart)
		optimize_init_argument(S, 1, nn)
		optimize_init_argument(S, 2, y)
		optimize_init_argument(S, 3, x1)
		optimize_init_which(S,"max")
		bhat1=optimize(S)
		
 		conv= optimize_result_converged(S)
        if (conv==0) {
		"warning: convergence not achieved, and now using the best RF estimates"
		}
		
		if (db==1){
			"reduced form parameters"
			bhat1
		}
		
		bhat1[1]=bhat1[1]+0.01	
		
		if (db==1){
			bo_logit(todo, bhat1, nn, y, x1, f, g, H)
			"initial parameters:"
			bstart
			
			"initial f"
			colsum(f)
			"initial g"
			colsum(g)
			delta=1.0d-5
			"numeric derivatives"
			for (i=1; i<=cols(bhat1); i++) {
				bb1=bhat1
				bb1[i]=bb1[i]+delta	
		     	bo_logit(todo, bb1, nn, y, x1, f1, g, H)
				
				colsum(f1-f)/delta	
			}
		}
				
		H1=optimize_result_Hessian(S)
		sco1=optimize_result_scores(S)	
		phat=Lambda(x1*bhat1')		
		
		bstart=J(1,cols(x2),0)
		if (ak>0) {
			kapstart=0.1
			bstart=bstart, (ln(kapstart)-ln(1-kapstart))
		}
					
		S=optimize_init()
		if (pd==0) {
			optimize_init_verbose(S,0)
			optimize_init_tracelevel(S, "none")              
		}
		optimize_init_evaluator(S, &bo_fct1v_logit())
		optimize_init_evaluatortype(S, "v1")
		optimize_init_technique(S,"nr")
		delt=J(1,cols(bstart),0.1)
		optimize_init_nmsimplexdeltas(S,delt)
		optimize_init_params(S,bstart)
		optimize_init_argument(S, 1, nn)
		optimize_init_argument(S, 2, y)
		optimize_init_argument(S, 3, x2)
		optimize_init_argument(S, 4, phat)
		optimize_init_argument(S, 5, ak)

		optimize_init_which(S,"max")
		bhat2=optimize(S)
		
		if (conv==0) {
		"warning: convergence not achieved, and now using the best coeff. estimates"
		}
		
		H2=optimize_result_Hessian(S)
		sco2=optimize_result_scores(S)

		bhat=bhat1,bhat2

        todo=1
        bo_fct1v_logit(todo, bhat2, nn, y, x2, phat, ak, f, g, H)
        df=colsum(g)
		                
        delta=1.0d-5
				
        H21=J(cols(bhat2),cols(bhat1),0)
        for (i=1; i<=cols(H1); i++) {
            bb1=bhat1
            bb1[i]=bb1[i]+delta
			phat=Lambda(x1*bb1')	
            bo_fct1v_logit(todo, bhat2, nn, y, x2, phat,ak, f, g, H)
            df1=colsum(g)
            H21[.,i]=(df1-df)':/delta
		}
		
		  
		H=(H1,J(rows(H1),cols(H2),0))\(H21,H2)
		sco=sco1,sco2
		
		H=-H
		if (pd==1) {
		sv=svdsv(H) 
		printf("\nRatio of largest to smallest singular value is %f \n",sv[1]/sv[cols(H)])
		}
		H=luinv(H)	
	
		if (sum(H,1)==.){
			printf("WARNING: H singular\n")
			printf("Call Bo at 609 216 5166\n\n")
			V=J(cols(sco),cols(sco),9.9999d99)
		}
		else {
			V=H*(sco'*sco)*H'
		}
		
		/* bhat=J(1,cols(xx)+1,1)*/
		vhat=V
		fopt=f
		
}

void bo_fct1v_logit(real scalar todo,
			real vector theta,
			real matrix nn,
			real matrix y,
			real matrix x,
			real matrix phat,
			real scalar ak,
			real vector f,
			real matrix g,
			real matrix H
)
/* calculates objective function ---- vector version  */
/* kappa is the probability that the right one is not in the set of matches*/
/* kappa=exp(theta)/(1+exp(theta)); hence theta=ln(kappa)-ln(1-kappa) */
{	
	real vector beta, xb, pr1, pr0, pp, fi,y0,y1 
	real scalar i, n1, n2, dfi, kappa1
	
	f=J(rows(nn),1,0)
	g=J(rows(nn),cols(theta),0)
	H=J(cols(theta),cols(theta),0)
	
	beta=theta
	kappa1=0
	if (ak>0) {
		beta=theta[1..cols(theta)-1]
		kappa1=Lambda(theta[cols(theta)])
	}
	
	xb=x*beta'
	pr1=Lambda(xb)
	pr0=1:-pr1
	y1=y:/phat
	y0=(1:-y):/(1:-phat)
	pp=1:/nn[.,3]

	for (i=1; i<=rows(nn); i++) {
		n1=nn[i,1]
		n2=nn[i,2]
		fi=(y1[n1..n2]'*pr1[n1..n2]+y0[n1..n2]'*pr0[n1..n2])*pp[i]
		f[i]=ln(fi*(1-kappa1)+kappa1)
		dfi=pp[i]:*(y1[n1..n2]-y0[n1..n2]):*pr1[n1..n2]:*pr0[n1..n2]
		if (ak>0) {
			g[i,.]=((1-kappa1):*dfi'x[n1..n2,.]:/(fi*(1-kappa1)+kappa1), (kappa1*(1-kappa1)) :*( -fi+1)  :/(fi*(1-kappa1)+kappa1))
		}
		else {
			g[i,.]=dfi'x[n1..n2,.]:/fi
		}			
	}
	if (todo >=2){
	H=0		
	}		
}


void bo_logit(real scalar todo,
			real vector theta,
			real matrix nn,
			real matrix y,
			real matrix x,
			real vector f,
			real matrix g,
			real matrix H
)
/* calculates logit objective function  */

{	
	real vector beta, xb, fi, pr1
	real scalar i, n1, n2
	real matrix dfi
	
	f=J(rows(nn),1,0)
	g=J(rows(nn),cols(theta),0)
	H=J(cols(theta),cols(theta),0)
	beta=theta
	xb=x*beta'
	pr1=Lambda(xb)
	
	fi=y:*log(pr1):+(1:-y):*log(1:-pr1)
	dfi=(y-pr1):*x 
/*	fi=log(y:*pr1:+(1:-y):*(1:-pr1)) */
/*	dfi = x:*(2*y:-1):*pr1:/(1:-y:+y:*exp(xb)) */

	for (i=1; i<=rows(nn); i++) {
		n1=nn[i,1]
		n2=nn[i,2]
		f[i]=colsum(fi[n1..n2])	
		g[i,.]=colsum(dfi[n1..n2,.])

	}

	if (todo>=2){
	H=0		
	}		
}

real vector function Lambda(real vector x)
{
	real vector ex
	ex=exp(x)
	return(ex:/(1:+ex))
}

void bo_OLS(real matrix y, real matrix x, real matrix bhat, real matrix S2, real vector yhat)
{
	real matrix xx, xy, e


	xy=x'*y
	xx=x'x
	xx=invsym(xx)
	bhat=xx*xy
	yhat=x*bhat
	e=y-yhat
	S2=(e'*e)/(rows(y)-1)
	bhat=bhat'
}

void bo_kill()
{
		"Deliberate kill"
		(1,2)*(1,2,3)
}

end


