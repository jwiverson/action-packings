SphericalFunctions:=function(G,H)
	#Input:		A group G and a subgroup H.
	#Output:	A table giving the values of the spherical functions for the action of
	#		G on G/H. Except for the first row and last column, the columns index double 
	#		cosets of H in G, and the rows index the spherical functions. An entry in 
	#		the interior of the table gives the value of a spherical function on a 
	#		double coset. The first row gives the sizes of the double cosets in the
	#		corresponding columns. The last column gives the degrees of the characters
	#		behind the spherical functions. The last entry of the first row is always 0.
	#Note:		The ordering of the rows agrees with the command
	#			ConstituentsOfCharacter(PermutationCharacter(G,H));
	#		Remember, however, that the first row is occupied by sizes. Hence, the
	#		i-th constituent of the permutation character corresponds to the (i+1)-st row
	#		of the table.

	local oh,tbl,perm,const,j,dcsts,reps,sizes,k,sigma,newReps,newSizes,positions,a,pos,h,table,chi,vals,val;

	oh:=Order(H);

	tbl:=CharacterTable(G);
	perm:=PermutationCharacter(G,H);
	const:=ShallowCopy(ConstituentsOfCharacter(perm));


	dcsts:=DoubleCosetRepsAndSizes(G,H,H);
	reps:=List(dcsts,x->x[1]);
	sizes:=List(dcsts,x->x[2]);

	#For each double coset representative, pre-compute the positions of the conjugacy classes 
	#for each element of the right coset of the inverse of the representative
	positions:=NewDictionary(reps[1],true);
	for a in reps do
		pos:=[];
		for h in H do
			Add(pos, PositionProperty( ConjugacyClasses( tbl ), C -> h*a^(-1) in C ));
		od;
		AddDictionary(positions,a,pos);
	od;

	#For each constituent and each double coset rep, use the conjugacy class possitions to compute 
	#the value of the spherical function
	table:=[Concatenation(sizes,[0])];
	for chi in const do
		vals:=[];
		for a in reps do
			pos:=LookupDictionary(positions,a); #positions of the elements of the right coset of a^-1
			Add(vals,Sum(pos,k->chi[k])/oh); #the value of the spherical function for chi at a
		od;
		Add(vals,DegreeOfCharacter(chi)); #the last column gives the degrees
		Add(table,vals);
	od;

	return table;
end;