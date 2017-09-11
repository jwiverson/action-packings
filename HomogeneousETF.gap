#	This file contains a companion method, HomogeneousETF, for Example 4.5 in the paper
#	"Optimal Line Packings From Finite Group Actions", by Joseph W. Iverson,
#	John Jasper, and Dustin G. Mixon. It produces the Gram matrix of any ETF
#	described in Table 1 or Table 2 of that paper.
#
#	WARNING:	This could take a long time, on the order of several hours.


ProjectiveReduction:=function(G)
	#Input:		A Gram matrix G of some frame
	#Output:	The Gram matrix of a projective reduction of the same frame

	local n,unique,i,row,fn,new,j,oldRow,phase;

	n:=Size(G);
	unique:=[];

	for i in [1..n] do
		#check to see if the i-th row is a multiple of one we've seen
		row:=G[i];;
		fn:=PositionProperty(row,x->x<>0); #position of first nonzero entry
	
		new:=true; #it will be false if it matches an old row
		for j in unique do
			oldRow:=G[j];
			if oldRow[fn] = 0 then
				#they definitely don't match
				continue;
			else
				phase:=row[fn]/oldRow[fn];
				if row = phase*oldRow then
					#we've seen this one
					new:=false;
					break;
				fi;
			fi;
		od;

		if new then
			AddSet(unique,i);
		fi;
	od;

	return G{unique}{unique};
end;




HomogeneousETF:=function(arg)
	#Input:		A list v=[d,n,m,t,k] indicating the desired dimension
	#		(d) and number of frame vectors (n), as well as the transitive
	#		permutation group TransitiveGroup(m,k) whose action on
	#		t points is supposed to be transitive and have some
	#		spherical functions corresponding to constituents of 
	#		multiplicity one, that produce a d x n ETF. A second,
	#		optional argument is a string, either "real" or "complex".
	#
	#Output:	The Gram matrix of an ETF as explained above. If the second
	#		argument is present, the output is guaranteed to have
	#		entries in the desired field.
	#
	#		This function returns fail if any of the following happen:
	#		1) TransitiveGroup(m,k) is not t-transitive.
	#		2) The list of spherical functions contains more than
	#		   20 candidates of multiplicity one.
	#		3) There is no ETF corresponding to the input data.

	local v,d,n,m,t,k,G,tbl,x0,i,H,chi,const,dcsts,reps,sizes,positions,a,table,deg,vals,pos,h,picky,success,it,sum,rank,absSqr,stabPos,stabSize,num,real,z,complex,C,D,c,A,x,R,y,p,dcst,j;

	v:=arg[1]; d:=v[1]; n:=v[2]; m:=v[3]; t:=v[4]; k:=v[5];

	#Create the permutation group and the stabilizer

	G:=TransitiveGroup(m,k);

	if Transitivity(G) < t then
		Print("Error: Group is not 2-transitive.\n");
		return fail;
	fi;

	x0:=[];
	for i in [1..t] do
		Add(x0,i);
	od;

	H:=Stabilizer(G,x0,OnTuples);


	#Create the permutation character, and collect constituents
	#of multiplicity one

	tbl:=CharacterTable(G);
	chi:=PermutationCharacter(G,H);
	const:=ConstituentsOfCharacter(chi);
	const:=Filtered(const, psi->ScalarProduct(chi,psi)=1);

	#Check to ensure there are no more than 20 candidates
	if Size(const) > 20 then
		Print("Error: Too many candidates.\n");
		return fail;
	fi;


	#Find double cosets and create the table of spherical functions
	dcsts:=DoubleCosetRepsAndSizes(G,H,H);
	reps:=List(dcsts,x->x[1]);
	sizes:=List(dcsts,x->x[2]);

	#For each double coset representative, pre-compute the positions of the
	#conjugacy classes for each element of the right coset of the inverse of 
	#the representative

	positions:=NewDictionary(reps[1],true);
	for a in reps do
		pos:=[];
		for h in H do
			Add(pos, PositionProperty( ConjugacyClasses( tbl ), C -> h*a^(-1) in C ));
		od;
		AddDictionary(positions,a,pos);
	od;

	#For each constituent and each double coset rep, use the conjugacy class
	#positions to compute the value of the spherical function

	table:=[sizes];
	for chi in const do
		deg:=DegreeOfCharacter(chi);
		vals:=[];
		for a in reps do
			pos:=LookupDictionary(positions,a); 
			#positions of the elements of the right coset of a^-1
			Add(vals,deg*Sum(pos,k->chi[k])/Size(H)); 
			#the value of the scaled spherical function for chi at a
		od;
		Add(table,vals);
	od;


	#Iterate through combinations of rows until a dxn ETF appears
	if Size(arg)>1 then
		picky:=true;
	else
		picky:=false;
	fi;

	success:=false;

	it:=IteratorOfCombinations([2..Size(table)]);

	for v in it do
		if Length(v) = 0 then
			#Ignore the empty set
			continue;
		fi;

		sum:=Sum(v,i->table[i]); 
		#the sum of the rows indexed by v

		rank:=sum[1];

		if rank <> d then
			#We have the wrong dimension
			continue;
		fi;

		absSqr:=List([1..Size(table[i])],j->sum[j]*ComplexConjugate(sum[j])); 
		#the squares of the absolute values of the inner products

		#Check to ensure we have the right size
		stabPos:=Positions(absSqr,rank^2); 
		#the positions corresponding to the projective stabilizer

		stabSize:=Sum(stabPos,j->table[1][j]); 
		#the size of the projective stabilizer

		num:=Size(G)/stabSize; 
		#the number of vectors in the frame corresponding to v

		if num <> n then
			#We have the wrong number of vectors
			continue;
		fi;

		#Check equiangularity
		if Size(AsSet(absSqr))=2 then
			#We found an ETF!
			if picky then
				#Check to see if its real or complex
				real:=true;
				for z in sum do
					if z <> ComplexConjugate(z) then
						#not real
						real:=false;
					fi;
				od;
				complex:=not(real);

				#Compare against what we wanted
				if arg[2]="real" and real then
					success:=true;
					break;
				fi;

				if arg[2]="complex" and complex then
					success:=true;
					break;
				fi;

			else
				#We don't care if it's real or not
				success:=true;
				break;
			fi;
		fi;
	od;


	#Report failure if we didn't find an ETF
	if not(success) then
		Print("Error: No such ETF.\n");
		return fail;
	fi;

	
	#Now we construct the projections for the rows indexed by v,
	#and return their sum


	#Find the right and double cosets for H in G
	C:=List(RightCosets(G, H), x -> Representative(x));
	num:=Size(C);
	D:=DoubleCosets(G,H,H);
	c:=Size(D);

	#Construct the adjacency matrices for the Schurian association scheme
	#Make empty adjacency matrices, and then fill them in
	A := []; #a list of adjacency matrices
	for i in [1..c] do
		Add(A,[]);
	od;

	for x in [1..num] do
		#Make an empty row for each adjacency matrix
		R := []; #a list of empty rows
		for i in [1..c] do
			Add(R,[]);
		od;

		#Fill in the rows with 0s and 1s
		for y in [1..num] do
			for i in [1..c] do
				if C[y]*C[x]^(-1) in D[i] then
					Add(R[i],1);
				else
					Add(R[i],0);
				fi;
			od;
		od;

		#Add the rows into the adjacency matrices
		for i in [1..c] do
			Add(A[i],R[i]);
		od;
	od;

	#Make the linear combination of adjacency matrices using the values in sum
	p:=ZeroMutable(A[1]);
	for i in [1..c] do
		#Find the position of the representative corresponding to D[i]
		for j in [1..Size(table[1])] do
			if reps[j] in D[i] then
				p:=p+(1/num)*sum[j]*A[i];
			fi;
		od;
	od;

	
	#Finally, projectively reduce the result
	return stabSize/Size(H)*ProjectiveReduction(p);
end;