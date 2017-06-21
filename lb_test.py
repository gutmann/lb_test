#!/usr/bin/env python
#
# ; % 2D Lattice Boltzmann (BGK) model of a fluid.
# ; %  c4  c3   c2  D2Q9 model. At each timestep, particle densities propagate
# ; %    \  |  /    outwards in the directions indicated in the figure. An
# ; %  c5 -c9 - c1  equivalent 'equilibrium' density is found, and the densities
# ; %    /  |  \    relax towards that state, in a proportion governed by omega.
# ; %  c6  c7   c8      Iain Haslam, March 2006.
# ;
# ; ;; constants
# function latb2d, F=F, bound=bound, slow=slow, ux=ux,uy=uy,  $
# 	maxsteps=maxsteps, minsteps=minsteps
#
# 	; set up constants for calculations
# 	omega=1.0d ; relaxation term
# 	density=1.0d
# 	t1=4.0/9 ; non-movement coefficient
# 	t2=1.0/9 ; moving perpendicular
# 	t3=1.0/36; moving diagonal
# 	c_squ=1.0/3
# 	nx=31 ; number of x grid cells
# 	ny=31 ; number of y grid cells
#
# 	if not keyword_set(maxsteps) then maxsteps = 4000
# 	if not keyword_set(minsteps) then minsteps = 100
#
# 	; F=repmat(density/9,[nx ny 9]); FEQ=F; msize=nx*ny; CI=[0:msize:msize*7];
# 	; F is the array of probability of particle movement in each direction
# 	; nx x ny x 9 possible directions of movement (or non-movement)
# 	if not keyword_set(F) then  $
# 		F = make_array(nx,ny,9, value=density/9)
# 	FEQ=f ; F equilibrium initially equals F
#
# 	;; matrix size
# 	msize=nx*ny
# 	; CI is an array of offsets into the matrix.  each element moves forward by one nx x ny grid
# 	CI=indgen(8)*msize
#
# 	; %a=repmat(-15:15,[31,1]);BOUND=(a.^2+a'.^2)<16;BOUND(1:nx,[1 ny])=1;
# 	; %BOUND=zeros(nx,ny);BOUND(1:nx,1)=1;%open channel
# 	; BOUND=rand(nx,ny)>0.7; %extremely porous random domain
# 	if not keyword_set(BOUND) then  $
# 		BOUND = (randomu(seed, nx,ny) gt 0.8)
#
# ;	BOUND=intarr(nx,ny)
# ;	bound[*,0]=1
# ;	bound[10,0:10]=1
# 	; ON=find(BOUND); %matrix offset of each Occupied Node
# 	ON=where(bound)
# 	; TO_REFLECT=[ON+CI(1) ON+CI(2) ON+CI(3) ON+CI(4) ...
# 	;             ON+CI(5) ON+CI(6) ON+CI(7) ON+CI(8)];
# 	TO_REFLECT=[ON+CI[0], ON+CI[1], ON+CI[2], ON+CI[3], $
# 		ON+CI[4], ON+CI[5], ON+CI[6], ON+CI[7]]
# 	; REFLECTED= [ON+CI(5) ON+CI(6) ON+CI(7) ON+CI(8) ...
# 	;             ON+CI(1) ON+CI(2) ON+CI(3) ON+CI(4)];
# 	REFLECTED= [ON+CI[4], ON+CI[5], ON+CI[6], ON+CI[7], $
# 		ON+CI[0], ON+CI[1], ON+CI[2], ON+CI[3]]
#
# 	; avu=1; prevavu=1; ts=0; deltaU=1e-7; numactivenodes=sum(sum(1-BOUND));
# 	avu=1
# 	prevavu=1
# 	ts=0
# 	deltaU=1e-7 *sqrt((nx/total(1-BOUND[0,*]))) ; increase inlet pressure more if there are fewer inlet nodes active
# 	numactivenodes=total(1-BOUND)
#
# 	shiftleft=[(indgen(nx-1)+1),0]
# 	shiftright=[nx-1,indgen(nx-1)]
# 	shiftup=[indgen(ny-1)+1,0]
# 	shiftdown=[ny-1,indgen(ny-1)]
# 	; while (ts<4000 & 1e-10<abs((prevavu-avu)/avu)) | ts<100
# 	while ((ts lt maxsteps) and (1e-10 lt abs((prevavu/avu)/avu))) or (ts lt minsteps) do begin
# 		;     % Propagate
# 		;     F(:,:,4)=F([2:nx 1],[ny 1:ny-1],4);F(:,:,3)=F(:,[ny 1:ny-1],3);
# 		F[*,*,3]=F[*,shiftdown,3]
# 		F[*,*,3]=F[shiftleft,*,3]
# 		F[*,*,2]=F[*,shiftdown,2]
# 		;     F(:,:,2)=F([nx 1:nx-1],[ny 1:ny-1],2);F(:,:,5)=F([2:nx 1],:,5);
# 		F[*,*,1]=F[*,shiftdown,1]
# 		F[*,*,1]=F[shiftright,*,1]
# 		F[*,*,4]=F[shiftleft,*,4]
# 		;     F(:,:,1)=F([nx 1:nx-1],:,1);F(:,:,6)=F([2:nx 1],[2:ny 1],6);
# 		F[*,*,0]=F[shiftright,*,0]
# 		F[*,*,5]=F[shiftleft,*,5]
# 		F[*,*,5]=F[*,shiftup,5]
# 		;     F(:,:,7)=F(:,[2:ny 1],7); F(:,:,8)=F([nx 1:nx-1],[2:ny 1],8);
# 		F[*,*,6]=F[*,shiftup,6]
# 		F[*,*,7]=F[shiftright,*,7]
# 		F[*,*,7]=F[*,shiftup,7]
# 		;     BOUNCEDBACK=F(TO_REFLECT); %Densities bouncing back at next timestep
# 		BOUNCEDBACK=F[TO_REFLECT]; %Densities bouncing back at next timestep
# 		;     DENSITY=sum(F,3);
# 		DENSITY=total(F,3)
# 		;     UX=(sum(F(:,:,[1 2 8]),3)-sum(F(:,:,[4 5 6]),3))./DENSITY;
# 		UX=(total(F[*,*,[0, 1, 7]],3)-total(F[*,*,[3,4,5]],3))/DENSITY
# 		;     UY=(sum(F(:,:,[2 3 4]),3)-sum(F(:,:,[6 7 8]),3))./DENSITY;
# 		UY=(total(F[*,*,[1,2,3]],3)-total(F[*,*,[5,6,7]],3))/DENSITY
# 		;     UX(1,1:ny)=UX(1,1:ny)+deltaU; %Increase inlet pressure
# 		UX[0,0:ny-1]+=deltaU; %Increase inlet pressure
# 		;     UX(ON)=0; UY(ON)=0; DENSITY(ON)=0;
# 		UX[ON]=0 & UY[ON]=0 & DENSITY[ON]=0
# 		;     U_SQU=UX.^2+UY.^2; U_C2=UX+UY; U_C4=-UX+UY; U_C6=-U_C2; U_C8=-U_C4;
# 		U_SQU=UX^2+UY^2 & U_C2=UX+UY & U_C4=-1*UX+UY & U_C6=-1*U_C2 & U_C8=-1*U_C4
# 		;     % Calculate equilibrium distribution: stationary
# 		;     FEQ(:,:,9)=t1*DENSITY.*(1-U_SQU/(2*c_squ));
# 		FEQ[*,*,8]=t1*DENSITY*(1-U_SQU/(2*c_squ))
# 		;     % nearest-neighbours
# 		;     FEQ(:,:,1)=t2*DENSITY.*(1+UX/c_squ+0.5*(UX/c_squ).^2-U_SQU/(2*c_squ));
# 		FEQ[*,*,0]=t2*DENSITY*(1+UX/c_squ+0.5*(UX/c_squ)^2-U_SQU/(2*c_squ))
# 		;     FEQ(:,:,3)=t2*DENSITY.*(1+UY/c_squ+0.5*(UY/c_squ).^2-U_SQU/(2*c_squ));
# 		FEQ[*,*,2]=t2*DENSITY*(1+UY/c_squ+0.5*(UY/c_squ)^2-U_SQU/(2*c_squ));
# 		;     FEQ(:,:,5)=t2*DENSITY.*(1-UX/c_squ+0.5*(UX/c_squ).^2-U_SQU/(2*c_squ));
# 		FEQ[*,*,4]=t2*DENSITY*(1-UX/c_squ+0.5*(UX/c_squ)^2-U_SQU/(2*c_squ));
# 		;     FEQ(:,:,7)=t2*DENSITY.*(1-UY/c_squ+0.5*(UY/c_squ).^2-U_SQU/(2*c_squ));
# 		FEQ[*,*,6]=t2*DENSITY*(1-UY/c_squ+0.5*(UY/c_squ)^2-U_SQU/(2*c_squ));
# 		;     % next-nearest neighbours
# 		;     FEQ(:,:,2)=t3*DENSITY.*(1+U_C2/c_squ+0.5*(U_C2/c_squ).^2-U_SQU/(2*c_squ));
# 		;     FEQ(:,:,4)=t3*DENSITY.*(1+U_C4/c_squ+0.5*(U_C4/c_squ).^2-U_SQU/(2*c_squ));
# 		;     FEQ(:,:,6)=t3*DENSITY.*(1+U_C6/c_squ+0.5*(U_C6/c_squ).^2-U_SQU/(2*c_squ));
# 		;     FEQ(:,:,8)=t3*DENSITY.*(1+U_C8/c_squ+0.5*(U_C8/c_squ).^2-U_SQU/(2*c_squ));
# 		FEQ[*,*,1]=t3*DENSITY*(1+U_C2/c_squ+0.5*(U_C2/c_squ)^2-U_SQU/(2*c_squ));
# 		FEQ[*,*,3]=t3*DENSITY*(1+U_C4/c_squ+0.5*(U_C4/c_squ)^2-U_SQU/(2*c_squ));
# 		FEQ[*,*,5]=t3*DENSITY*(1+U_C6/c_squ+0.5*(U_C6/c_squ)^2-U_SQU/(2*c_squ));
# 		FEQ[*,*,7]=t3*DENSITY*(1+U_C8/c_squ+0.5*(U_C8/c_squ)^2-U_SQU/(2*c_squ));
# 		;     F=omega*FEQ+(1-omega)*F;
# 		F=omega*FEQ+(1-omega)*F
# 		;     F(REFLECTED)=BOUNCEDBACK;
# 		F[REFLECTED]=BOUNCEDBACK
# 		;     prevavu=avu;avu=sum(sum(UX))/numactivenodes; ts=ts+1;
# 		prevavu=avu
# 		avu=total(UX)/numactivenodes
# 		ts++
# 		if (ts mod 10) eq 0 and keyword_set(slow) then begin
# 			velovect,ux[1:*,*],uy[1:*,*],length=1
# ;			wait, 0.1
# 		endif
# 	endwhile
# 	; end
# 	; figure;colormap(gray(2));image(2-BOUND');hold on;
# 	; quiver(2:nx,1:ny,UX(2:nx,:)',UY(2:nx,:)');
# 	; title(['Flow field after ',num2str(ts),'\deltat']);xlabel('x');ylabel('y');
# ;	window, xs=330,ys=330
# ;	for i=0,8 do $
# ;		tvscl, congrid(F[*,*,i],100,100),i
# ;	velovect,ux[1:*,*],uy[1:*,*],length=1
# 	;; see also vel and plot_field
# 	return, F
# end
#
# pro latmovie, nops=nops, b=b
#
# 	if not keyword_set(nops) then old=setupplot(filename='latbmovie.ps',xs=8,ys=7)
# 	!p.multi=[0,1,1]
# 	usersym, [1,1,-1,-1],[-1,1,1,-1],/fill
# 	density=1.0d
# 	nx=31
# 	ny=31
# 	ss=1.08
# 	if keyword_set(nops) then ss=3
#
# 	F = make_array(nx,ny,9, value=density/9)
# 	BOUND = (randomu(seed, nx,ny) gt 0.8)
# 	BOUND=fltarr(nx,ny)
# 	BOUND[0,5:*]=1
# 	BOUND[*,0]=1
# 	movingoffset=nx/2
# 	BOUND[movingoffset,4:*]=1
#
# 	;; inital equilibrium LB calc
# 	F=latb2d(F=F,bound=bound, ux=ux,uy=uy, maxsteps=2000)
#
# 	velovect, ux[1:*,*],uy[1:*,*],length=1
# 	tmp=where(bound)
# 	oplot, (tmp mod nx)-1, tmp / nx, psym=8,symsize=ss
#
# 	for i=0,100 do begin
# 		bot=i mod ny
# 		top=(bot+3) mod ny
# 		if top gt bot then begin
# 			bound[movingoffset,*]=1
# 			bound[movingoffset,bot:top]=0
# 		endif else begin
# 			bound[movingoffset,*]=0
# 			bound[movingoffset,top:bot]=1
# 		endelse
# 		F=latb2d(F=F,bound=bound,ux=ux,uy=uy, maxsteps=200, minsteps=100)
# 		tmp=where(bound)
# 		if not keyword_set(b) then begin
# 			velovect, ux[1:*,*],uy[1:*,*],length=1.5
# 			oplot, (tmp mod nx)-1, tmp / nx, psym=8,symsize=ss
# 		endif else begin
# 			vel, ux[1:*,*],uy[1:*,*],length=0.4
# 			oplot, ((tmp mod nx))/float(nx), ((tmp/nx)+0.5)/float(ny), psym=8,symsize=ss
# 		endelse
# ;		plot_field, ux[1:*,*],uy[1:*,*],length=1
# ;		oplot, ((tmp mod nx)-1)/float(nx), ((1+tmp) / nx)/float(ny), psym=8,symsize=ss
# 		if !d.name eq 'X' then wait, 0.1
# 	endfor
# 	if not keyword_set(nops) then resetplot, old
# end

def main():
    print("This file contains older IDL and Matlab Lattice Boltzmann code. ")

if __name__ == '__main__':
    main()
