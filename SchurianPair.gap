#This file is a slight modification of the method IsSchurian, created by Prof. Akihide Hanaki.
#At the time of writing, that code was available as part of a much larger collection of 
#functions for analyzing association schemes in GAP, at 
#  http://math.shinshu-u.ac.jp/~hanaki/as/
#It relies, in turn, on the graph automorphism methods of the GRAPE package, which ultimately 
#uses nauty to do the heavy lifting.

LoadPackage("grape);

SchurianPair := function(R)
	#Input:		An association scheme R as stored by Hanaki's method: if A[i] is the
	#		i-th adjacency matrix, then R:=Sum([0..d],i->i*A[i]);
	#Output:	If the scheme is Schurian, this function returns a list [G,H],
	#		consisting a group G and a subgroup H, such that Hanaki's function
	#		SchurianScheme(G,H) returns an association scheme isomorphic
	#		to R. If the scheme is not Schurian, this function returns "fail"

	local G, adj, gp, gr, n, x, y, i;
    
	adj := AdjacencyMatrices(R);
	n := Length(R);
	G := SymmetricGroup(n);
	for i in [2..(Length(adj) - 1)] do
		gr := Graph(Group(()), [1..n], OnPoints, function(x,y) return adj[i][x][y]=1; end);
		gp := AutomorphismGroup(gr);
		G := Intersection(G, gp);
		if not(IsTransitive(G, [1..n])) then
			return fail;
		fi;
	od;
  
	if not( (Length(Orbits(Stabilizer(G, 1), [1..n])) = Length(adj)) ) then
		return fail;
	else
		return [G,Stabilizer(G,1)];
	fi;
end;