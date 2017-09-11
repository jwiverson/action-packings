SphericalProjections:=function(G,H)
	#Input:		A group G and a subgroup H.
	#Output:	A list containing the projections in the Schurian association scheme associated
	#		with spherical functions for the action of G on G/H.
	#Note:		The output is compatible with the function SchurianScheme in the file of
	#		association scheme commands produced by Akihide Hanaki, available at
	#			http://math.shinshu-u.ac.jp/~hanaki/as/gap/association_scheme.gap
	#		In other words, the projection we produce here lies in the scheme produced by that function.
	#		In fact, we modified some of the code from that file to produce the association scheme
	#		as part of this function.
	#Note 2:	The ordering of the projections agrees with the command
	#			ConstituentsOfCharacter(PermutationCharacter(G,H));
	#		In other words, the i-th constituent builds the i-th projection in our list.


	#Build the association scheme
	local C,n,D,d,A,i,x,y,R,tbl,perm,const,positions,dcst,pos,e,chi,m,c;

	C := List(RightCosets(G, H), x -> Representative(x));
	n := Size(C);
	D := DoubleCosets(G, H, H);
	d := Size(D);

	#Make empty adjacency matrices, and then fill them in
	A := []; #a list of adjacency matrices
	for i in [1..d] do
		Add(A,[]);
	od;

	for x in [1..n] do
		#Make an empty row for each adjacency matrix
		R := []; #a list of empty rows
		for i in [1..d] do
			Add(R,[]);
		od;

		#Fill in the rows with 0s and 1s
		for y in [1..n] do
			for i in [1..d] do
				if C[y]*C[x]^(-1) in D[i] then
					Add(R[i],1);
				else
					Add(R[i],0);
				fi;
			od;
		od;

		#Add the rows into the adjacency matrices
		for i in [1..d] do
			Add(A[i],R[i]);
		od;
	od;


	#Make the character table of G and compute the constituents of the permutation character
	tbl:=CharacterTable(G);
	perm:=PermutationCharacter(G,H);
	const:=ConstituentsOfCharacter(perm);

	#Pre-compute the conjugacy class positions of every double coset, to speed computation
	positions:=[];
	for dcst in D do
		pos:=[];
		for x in dcst do
			Add(pos, PositionProperty( ConjugacyClasses(tbl), C -> x in C ) );
		od;
		Add(positions,pos);
	od;

	#Compute the projections
	e:=[];
	for chi in const do
		m:=DegreeOfCharacter(chi);
		c:=[]; #a list of coefficients
		for i in [1..d] do
			Add(c, m/n/Size(D[i])*Sum(positions[i],x->chi[x]) );
		od;
		Add(e, Sum([1..d],i->c[i]*A[i]) );
	od;

	return e;
end;
