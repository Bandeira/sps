#include "spectral_alignment.h"
#include "aminoacid.h"
#include <fstream>

#include <vector>
#include <algorithm>
#include <cmath>

namespace SPECNETS_EXTRA {
	const double infinity = 1e10;
	const int None = -1;
	
	const double mass_tolerance = 0.5;
	
	using std::vector;
	using std::pair;
	using std::max;
	using std::min;
	using std::max_element;
	
	enum celltype {cell_prefix1, cell_suffix1, cell_prefix2, cell_suffix2,
		       cell_prefix1_L, cell_suffix1_R, cell_prefix2_L, cell_suffix2_R,
		       cell_D1, cell_D2, cell_D3,
		       cell_M1_L, cell_M1_R, cell_M2_L, cell_M2_R,
		       cell_M3_L, cell_M3_R, INVALID};
	
	typedef vector<vector<double> > vector2;
	typedef vector<vector<vector<double> > > vector3;
	typedef pair<vector<int>, vector<int> > pair2;
	
	
	//////////////////////////////////////////////////////////////////////////////
	pair<double, pair2> align_simple(double parent_mass1, double parent_mass2,
					 vector<double> peaks,
					 vector<double> peaks2,
					 vector<int> & common,
					 vector<int> & common2,
					 vector<double> & common_scores,
					 vector<double> & common2_scores,
					 vector<int> & prev,
					 vector<int> & next,
					 vector<vector<int> > & left_jumps,
					 vector<vector<int> > & right_jumps,
					 double ptm_penalty,
					 double ambiguous_penalty)
	// parent_mass1 : parent mass of first spectrum
	// parent mass2 : parent mass of 2nd spectrum.
	// Note: parent_mass2 can be either bigger or smaller than parent_mass1
	//
	// peaks        : All the peaks m in the 1st spectrum such that
	//                either the peak m or the peak m+(parent_mass2-parent_mass1)
	// peaks2       : All the peaks of the 2st spectrum
	//
	// common       : common[i] is the index j of the peak in spectrum2 whose mass
	//                is equal to peaks[i].
	//                If no such peak exists, common[i] = None
	// common2      : common[i] is the index j of the peak in spectrum2 whose mass
	//                is equal to peaks[i]+(parent_mass2-parent_mass1).
	// common_scores: common_scores[i] is the score of the vertex
	//                (peaks[i],peaks[i]) in the alignment graph.
	//                If the vertex doesn't exist, then common_scores[i]=-infinity
	// common2_scores: common_scores2[i] is the score of the vertex
	//                (peaks[i],peaks[i]+(parent_mass2-parent_mass1).
	//
	// Let  T = the mass tolerance (e.g. 0.5)
	//
	// Let J be a set of masses (possibly empty), and M be a mass.
	// We say that a mass jump m is valid if either |m-m'| < T for some m' from J,
	// or m >= M.
	// Assume that if J is not empty, then min(J) >= 57 and max(J) < M
	//
	// Define
	//   delta2 = M-T
	//   
	// prev         : prev[i] is the maximum index j such that
	//                peaks[j] < peaks[i]-delta2
	//                If no such index exists, prev[i] = None
	// next         : next[i] is the minimum index j such that
	//                peaks[j] > peaks[i]+delta2
	//
	// left_jumps   : left_jumps[i] is a list of all the indices j such that
	//   (1) peaks[i]-delta2 <= peaks[j] < peaks[i]-(57-T)
	//   (2) peaks[i]-peaks[j] is a valid mass jump
	//   If J is empty, then left_jumps[i] is an empty list
	// right_jumps   : right_jumps[i] is a list of all the indices j such that
	//   (1) peaks[i]+57-T < peaks[j] <= peaks[i]+delta2
	//   (2) peaks[j]-peaks[i] is a valid mass jump
	//   If J is empty, then right_jumps[i] is an empty list
	//
	// ptm_penalty:   penalty for having a PTM (negative number)
	// ambiguous_penalty: penalty if the first peak i before (or after) the
	//                    modification satisfies
	//                        common[i] != None and common2[i] != None
	{
	    int n = common.size();
	
	    vector<double> prefix(n, -infinity);
	    vector<double> suffix(n, -infinity);
	    vector<double> prefix_L(n, -infinity);
	    vector<double> suffix_R(n, -infinity);
	
	    // compute prefix
	    for (int i = 0; i < n; ++i) {
	        if (common[i] != None) {
	            double prev_score = 0.0;
		    for (vector<int>::const_iterator it = left_jumps[i].begin();
			 it != left_jumps[i].end(); ++it) {
			prev_score = max(prev_score, prefix[*it]);
		    }
	            int i2 = prev[i];
	            if (i2 != None) {
	                prev_score = max(prev_score, prefix_L[i2]);
		    }
	            prefix[i] = common_scores[i] + prev_score;
		}
	
		if (i != 0)
	            prefix_L[i] = max(prefix_L[i-1], prefix[i]);
	        else
	            prefix_L[i] = prefix[i];
	    }
	
	    // compute suffix
	    for (int j = n-1; j >= 0; --j) {
	        if (common2[j] != None) {
	            double next_score = 0.0;
		    for (vector<int>::const_iterator it = right_jumps[j].begin();
			 it != right_jumps[j].end(); ++it) {
			next_score = max(next_score, suffix[*it]);
		    }
	            int j2 = next[j];
	            if (j2 != None) {
	                next_score = max(next_score, suffix_R[j2]);
		    }
		    suffix[j] = common2_scores[j] + next_score;
		}
		if (j != n-1)
	            suffix_R[j] = max(suffix_R[j+1], suffix[j]);
	        else
	            suffix_R[j] = suffix[j];
	    }
	
	    double best_score = 0.0;
	    int best_i = None;
	    int best_j = None;
	    celltype best_t = INVALID;
	    for (int i = 0; i < n; ++i) {
		double penalty = ptm_penalty;
		if (common2[i] != None)
		    penalty += ambiguous_penalty;
		for (vector<int>::const_iterator it = right_jumps[i].begin();
		     it != right_jumps[i].end(); ++it) {
		    int j = *it;
		    double penalty2 = penalty;
		    if (common[j] != None)
			penalty2 += ambiguous_penalty;
		    if (best_score < prefix[i]+suffix[j]+penalty2) {
			best_i = i;
	                best_j = j;
	                best_t = cell_suffix2;
	                best_score = prefix[i]+suffix[j]+penalty2;
		    }
		}
		int j = next[i];
		if (j != None) {
		    if (best_score < prefix[i]+suffix_R[j]+penalty) {
			best_i = i;
			best_j = j;
			best_t = cell_suffix2_R;
			best_score = prefix[i]+suffix_R[j]+penalty;
		    }
		}
	    }
	
	    for (int i = 0; i < n; ++i) {
		double penalty = 0;
		if (common2[i] != None)
		    penalty += ambiguous_penalty;
		if (best_score < prefix[i]+penalty) {
	            best_i = i;
		    best_j = None;
		    best_t = cell_suffix2;
	            best_score = prefix[i]+penalty;
		}
	    }
	    for (int j = 0; j < n; ++j) {
		double penalty = 0;
		if (common[j] != None)
		    penalty += ambiguous_penalty;
	        if (best_score < suffix[j]+penalty) {
	            best_i = None;
	            best_j = j;
	            best_t = cell_suffix2;
	            best_score = suffix[j]+penalty;
		}
	    }
	
	    vector<int> path1;
	    vector<int> path2;
	
	    if (best_t == INVALID)
		return pair<double, pair2>(best_score/parent_mass1, pair2(path1, path2));
	
	    int i = best_i;
	    while (i != None) {
		path1.insert(path1.begin(), i);
	        double prev_score = 0.0;
	        int index = 0;
	        int next_i = 0;
		for (vector<int>::const_iterator it = left_jumps[i].begin();
		     it != left_jumps[i].end(); ++it) {
		    int i2 = *it;
	            if (prefix[i2] > prev_score) {
	                prev_score = prefix[i2];
	                index = 1;
	                next_i = i2;
		    }
		}
	        int i2 = prev[i];
	        if (i2 != None && prefix_L[i2] > prev_score) {
	            prev_score = prefix_L[i2];
	            index = 2;
	            next_i = i2;
		}
	
	        i = next_i;
		if (index == 0) {
	            break;
	        } else if (index == 1) {
	            // do nothing
		} else {
	            while (prefix_L[i] != prefix[i])
	                i -= 1;
		}
	    }
	
	    int j = best_j;
	    if (j != None && best_t == cell_suffix2_R) {
	        while (suffix_R[j] != suffix[j]) {
	            j += 1;
		}
	    }
	
	    while (j != None) {
		path2.push_back(j);
	        double next_score = 0.0;
	        int index = 0;
	        int next_j = 0;
		for (vector<int>::const_iterator it = right_jumps[j].begin();
		     it != right_jumps[j].end(); ++it) {
		    int j2 = *it;
	            if (suffix[j2] > next_score) {
	                next_score = suffix[j2];
	                index = 1;
	                next_j = j2;
		    }
		}
	        int j2 = next[j];
	        if (j2 != None && suffix_R[j2] > next_score) {
	            next_score = suffix_R[j2];
	            index = 2;
	            next_j = j2;
		}
	
	        j = next_j;
	        if (index == 0) {
	            break;
		} else if (index == 1) {
	            // do nothing
		} else {
	            while (suffix_R[j] != suffix[j]) {
	                j += 1;
		    }
		}
	    }
	
	    return pair<double, pair2>(best_score/parent_mass1, pair2(path1, path2));
	
	}
	
	//////////////////////////////////////////////////////////////////////////////
	
	pair<double, pair2> find_paths(double parent_mass1, double parent_mass2,
				       vector<double> & common_scores,
				       vector<double> & common2_scores,
				       double ptm_penalty,
				       vector<int> & prev, vector<int> & next,
				       vector<int> & prev2, vector<int> & next2,
				       vector<vector<int> > & left_jumps,
				       vector<vector<int> > & right_jumps,
				       vector<vector<int> > & left_jumps2,
				       vector<vector<int> > & right_jumps2,
				       vector<vector<int> > & left_neighbors,
				       vector<vector<int> > & right_neighbors,
				       int best_i, int best_j, int best_s,
				       celltype best_t, double best_score,
				       vector2 & D1, vector3 & D2, vector3 & D3,
				       vector2 & D2max, vector2 & D3max,
				       vector2 & M1_L, vector2 & M1_R,
				       vector2 & M2_L, vector3 & M2_R,
				       vector3 & M3_L, vector2 & M3_R,
				       vector<double> & prefix1, vector2 & suffix1,
				       vector2 & prefix2, vector<double> & suffix2,
				       vector<double> & prefix1_L,
				       vector<double> & suffix1max, vector<double> & suffix1_R,
				       vector<double> & prefix2max, vector<double> & prefix2_L,
				       vector<double> & suffix2_R);
	
	//////////////////////////////////////////////////////////////////////////////
	
	int max_ind(vector<double> & v)
	{
	    return max_element(v.begin(), v.end())-v.begin();
	}
	
	//////////////////////////////////////////////////////////////////////////////
	
	pair<double, pair2> align(double parent_mass1, double parent_mass2,
				  vector<double> peaks,
				  vector<double> peaks2,
				  vector<int> & common,
				  vector<int> & common2,
				  vector<double> & common_scores,
				  vector<double> & common2_scores,
				  vector<int> & prev,
				  vector<int> & next,
				  vector<int> & prev2,
				  vector<int> & next2,
				  vector<vector<int> > & left_jumps,
				  vector<vector<int> > & right_jumps,
				  vector<vector<int> > & left_jumps2,
				  vector<vector<int> > & right_jumps2,
				  vector<vector<int> > & left_neighbors,
				  vector<vector<int> > & right_neighbors,
				  double same_vertex_penalty,
				  double ptm_penalty)
	
	// parent_mass1 : parent mass of first spectrum
	// parent mass2 : parent mass of 2nd spectrum. We assume that
	//                parent_mass2 >= parent_mass1
	// peaks        : All the peaks m in the 1st spectrum such that
	//                either the peak m or the peak m+(parent_mass2-parent_mass1)
	// Note: peaks must be symmetric, namely, peaks[i]+peaks[n-1-i]=parent_mass1+18
	//       for all i, where n is the length of peaks
	//
	// peaks2       : All the peaks of the 2st spectrum
	// common       : common[i] is the index j of the peak in spectrum2 whose mass
	//                is equal to peaks[i].
	//                If no such peak exists, common[i] = None
	// common2      : common[i] is the index j of the peak in spectrum2 whose mass
	//                is equal to peaks[i]+(parent_mass2-parent_mass1).
	// common_scores: common_scores[i] is the score of the vertex
	//                (peaks[i],peaks[i]) in the alignment graph.
	//                If the vertex doesn't exist, then common_scores[i]=-infinity
	// common2_scores: common_scores2[i] is the score of the vertex
	//                (peaks[i],peaks[i]+(parent_mass2-parent_mass1).
	//
	//
	// Let  T = the mass tolerance (e.g. 0.5)
	//
	// Let J be a set of masses (possibly empty), and M be a mass.
	// We say that a mass jump m is valid if either |m-m'| < T for some m' from J,
	// or m >= M.
	// Assume that if J is not empty, then min(J) >= 57 and max(J) < M
	//
	// Define
	//   delta0 = parent_mass2-parent_mass1
	//   delta = min(max(delta0+T, 57-T),114-T)
	//   delta2 = M-T
	//   
	// prev         : prev[i] is the maximum index j such that
	//                peaks[j] < peaks[i]-delta2
	//                If no such index exists, prev[i] = None
	// next         : next[i] is the minimum index j such that
	//                peaks[j] > peaks[i]+delta2
	// prev2        : prev2[i] is the maximum index j such that
	//                peaks[j] < peaks[i]-max(delta, delta2).
	// next2        : next[i] is the minimum index j such that
	//                peaks[j] > peaks[i]+max(delta, delta2).
	//
	// left_jumps   : left_jumps[i] is a list of all the indices j such that
	//   (1) peaks[i]-delta2 <= peaks[j] < peaks[i]-(57-T)
	//   (2) peaks[i]-peaks[j] is a valid mass jump
	//   If J is empty, then left_jumps[i] is an empty list
	// right_jumps   : right_jumps[i] is a list of all the indices j such that
	//   (1) peaks[i]+57-T < peaks[j] <= peaks[i]+delta2
	//   (2) peaks[j]-peaks[i] is a valid mass jump
	//   If J is empty, then right_jumps[i] is an empty list
	//
	// left_jumps2  : left_jumps2[i] is a list of all the indices j such that
	//   (1) peaks[i]-delta2 <= peaks[j] <= peaks[i]-delta
	//   (2) peaks[i]-peaks[j] is a valid mass jump
	//   If delta >= delta2 then left_jumps2[i] is an empty list (follows from (1))
	// right_jumps2  : right_jumps2[i] is a list of all the indices j such that
	//   (1) peaks[i]+delta < peaks[j] <= peaks[i]+delta2
	//   (2) peaks[j]-peaks[i] is a valid mass jump
	//   If delta >= delta2 then right_jumps2[i] is an empty list (follows from (1))
	//
	// left_neighbors: left_neighbors[i] is a list of all indices j such that
	//   (1) peaks[i]-delta <= peaks[j] < peaks[i]-(57-T)
	//   (2) peaks[i]-peaks[j] is a valid mass jump
	// right_neighbors: right_neighbors[i] is a list of all indices j such that
	//   (1) peaks[i]+57-T < peaks[j] <= peaks[i]+delta
	//   (2) peaks[j]-peaks[i] is a valid mass jump
	//
	// same_vertex_penalty: penalty for using the same vertex twice 
	//                      (negative number)
	// ptm_penalty:   penalty for having a PTM (negative number)
	{
	    int n0 = common.size();
	    int N = n0-1;
	    int n = (n0+1)/2;
	
	    vector<double> prefix1(n, -infinity);
	    vector2 suffix1(n);
	    vector2 prefix2(n);
	    vector<double> suffix2(n, -infinity);
	
	    vector<double> prefix1_L(n, -infinity);
	    vector<double> suffix1max(n, -infinity);
	    vector<double> suffix1_R(n, -infinity);
	    vector<double> prefix2max(n, -infinity);
	    vector<double> prefix2_L(n, -infinity);
	    vector<double> suffix2_R(n, -infinity);
	
	    vector2 D1(n, vector<double>(n, -infinity));
	    vector3 D2(n, vector2(n));
	    vector3 D3(n, vector2(n));
	    vector2 D2max(n, vector<double>(n, -infinity));
	    vector2 D3max(n, vector<double>(n, -infinity));
	    vector2 M1_L(n, vector<double>(n, -infinity));
	    vector2 M1_R(n, vector<double>(n, -infinity));
	    vector2 M2_L(n, vector<double>(n, -infinity));
	    vector3 M2_R(n, vector2(n));
	    vector3 M3_L(n, vector2(n));
	    vector2 M3_R(n, vector<double>(n, -infinity));
	
	    // compute prefix1
	    for (int i = 0; i < n; ++i) {
	        if (common[i] != None) {
	            double prev_score = 0.0;
		    for (vector<int>::const_iterator it = left_jumps[i].begin();
			 it != left_jumps[i].end(); ++it) {
			prev_score = max(prev_score, prefix1[*it]);
		    }
	            int i2 = prev[i];
	            if (i2 != None) {
	                prev_score = max(prev_score, prefix1_L[i2]);
		    }
	            prefix1[i] = common_scores[i] + prev_score;
		}
	
		if (i != 0)
	            prefix1_L[i] = max(prefix1_L[i-1], prefix1[i]);
	        else
	            prefix1_L[i] = prefix1[i];
	    }
	
	    // compute suffix1
	    for (int j = 0; j < n; ++j) {
	        int l = right_neighbors[N-j].size();
	        suffix1[j].resize(l+1, -infinity);
	
	        if (common[N-j] != None) {
	            // s = 0
	            double next_score = 0.0;
		    for (vector<int>::const_iterator it = right_jumps[N-j].begin();
			 it != right_jumps[N-j].end(); ++it) {
			next_score = max(next_score, suffix2[N-*it]+ptm_penalty);
		    }
		    int j2 = next[N-j];
		    if (j2 != None) {
			next_score = max(next_score, suffix2_R[N-j2]+ptm_penalty);
		    }
		    for (vector<int>::const_iterator it = right_jumps2[N-j].begin();
			 it != right_jumps2[N-j].end(); ++it) {
	                next_score = max(next_score, suffix1max[N-*it]);
		    }
	            j2 = next2[N-j];
	            if (j2 != None) {
	                next_score = max(next_score, suffix1_R[N-j2]);
		    }
	
	            suffix1[j][0] = common_scores[N-j] + next_score;
		    // s > 0
		    for (int s = 0; s < l; ++s) {
			int j2 = right_neighbors[N-j][s];
	                if (common[j2] == None)
	                    continue;
			double penalty = 0.0;
			double y = peaks2[common[N-j]]+peaks2[common[j2]];
			if (fabs(y-(parent_mass2+18)) < mass_tolerance)
			    penalty = same_vertex_penalty;
	                double next_score = suffix1max[N-j2];
	                suffix1[j][s+1] = common_scores[N-j] + next_score + penalty;
		    }
		}
	
	        suffix1max[j] = *max_element(suffix1[j].begin(), suffix1[j].end());
	        if (j != 0)
	            suffix1_R[j] = max(suffix1_R[j-1], suffix1max[j]);
	        else
	            suffix1_R[j] = suffix1max[j];
	    }
	
	    // compute prefix2
	    for (int i = 0; i < n; ++i) {
		int l = left_neighbors[i].size();
		prefix2[i].resize(l+1 ,-infinity);
	
	        if (common2[i] != None) {
		    // s = 0
		    double prev_score = 0.0;
		    for (vector<int>::const_iterator it = left_jumps[i].begin();
			 it != left_jumps[i].end(); ++it) {
			prev_score = max(prev_score, prefix1[*it]+ptm_penalty);
		    }
		    int i2 = prev[i];
		    if (i2 != None)
			prev_score = max(prev_score, prefix1_L[i2]+ptm_penalty);
		    for (vector<int>::const_iterator it = left_jumps2[i].begin();
			 it != left_jumps2[i].end(); ++it) {
	                prev_score = max(prev_score, prefix2max[*it]);
		    }
	            i2 = prev2[i];
	            if (i2 != None) {
	                prev_score = max(prev_score, prefix2_L[i2]);
		    }
	            prefix2[i][0] = common2_scores[i] + prev_score;
		    // s > 0
		    for (int s = 0; s < l; ++s) {
			int i2 = left_neighbors[i][s];
	                if (common2[i2] == None)
	                    continue;
			double penalty = 0.0;
			double y = peaks2[common2[i]]+peaks2[common2[i2]];
			if (fabs(y-(parent_mass2+18)) < mass_tolerance)
			    penalty = same_vertex_penalty;
	                prev_score = prefix2max[i2];
	                prefix2[i][s+1] = common2_scores[i] + prev_score + penalty;
		    }
		}
	
	        prefix2max[i] = *max_element(prefix2[i].begin(), prefix2[i].end());
	        if (i != 0)
	            prefix2_L[i] = max(prefix2_L[i-1], prefix2max[i]);
	        else
	            prefix2_L[i] = prefix2max[i];
	    }
	
	    // compute suffix2
	    for (int j = 0; j < n; ++j) {
	        if (common2[N-j] != None) {
	            double next_score = 0.0;
		    for (vector<int>::const_iterator it = right_jumps[N-j].begin();
			 it != right_jumps[N-j].end(); ++it) {
			next_score = max(next_score, suffix2[N-*it]);
		    }
	            int j2 = next[N-j];
	            if (j2 != None) {
	                next_score = max(next_score, suffix2_R[N-j2]);
		    }
		    suffix2[j] = common2_scores[N-j] + next_score;
		}
		if (j != 0)
	            suffix2_R[j] = max(suffix2_R[j-1], suffix2[j]);
	        else
	            suffix2_R[j] = suffix2[j];
	    }
	
	    // compute D1/M1_R
	    for (int i = 0; i < n; ++i) {
	        if (common[i] == None) {
	            if (i != 0)
	                M1_L[i] = M1_L[i-1];
	            continue;
		}
	
		for (int j = 0; j < n; ++j) {
	            if (common2[N-j] == None) {
	                if (i != 0)
	                    M1_L[i][j] = M1_L[i-1][j];
	                if (j != 0)
	                    M1_R[i][j] = M1_R[i][j-1];
	                continue;
		    }
	
	            double penalty = 0.0;
	            double x = peaks[i]+peaks[N-j];
		    double y = peaks2[common[i]]+peaks2[common2[N-j]];
		    if (fabs(x-(parent_mass1+18)) < mass_tolerance ||
			fabs(y-(parent_mass2+18)) < mass_tolerance) 
			penalty = same_vertex_penalty;
	
		    if (i >= j) {
			double prev_score = suffix2[j]+ptm_penalty;
			// If we choose suffix2[j], then we pay the PTM penalty -
	                //  the penalty is paid only once!
			for (vector<int>::const_iterator it = left_jumps[i].begin();
			     it != left_jumps[i].end(); ++it) {
			    prev_score = max(prev_score, D1[*it][j]);
			}
	                int i2 = prev[i];
	                if (i2 != None)
	                    prev_score = max(prev_score, M1_L[i2][j]);
	                D1[i][j] = common_scores[i] + prev_score + penalty;
		    } else {
	                double next_score = prefix1[i]+ptm_penalty;
			for (vector<int>::const_iterator it = right_jumps[N-j].begin();
			     it != right_jumps[N-j].end(); ++it) {
			    next_score = max(next_score, D1[i][N-*it]);
			}
	                int j2 = next[N-j];
	                if (j2 != None) {
	                    next_score = max(next_score, M1_R[i][N-j2]);
			}
	                D1[i][j] = common2_scores[N-j] + next_score + penalty;
		    }
	
		    // Compute M1_L
	            if (i != 0)
	                M1_L[i][j] = max(M1_L[i-1][j], D1[i][j]);
	            else
	                M1_L[i][j] = D1[i][j];
	
		    // Compute M1_R
	            if (j != 0)
	                M1_R[i][j] = max(M1_R[i][j-1], D1[i][j]);
	            else
	                M1_R[i][j] = D1[i][j];
		}
	    }
	
	    // compute D2/D2max/M2_L
	    for (int i = 0; i < n; ++i) {
		int l = left_neighbors[i].size();
		for (int j = 0; j < n; ++j) {
	            D2[i][j].resize(l+1, -infinity);
	            M2_R[i][j].resize(l+1, -infinity);
		}
		if (common2[i] == None) {
	            if (i > 0)
	                M2_L[i] = M2_L[i-1];
	            continue;
		}
	
		for (int j = 0; j < n; ++j) {
	            if (common2[N-j] == None) {
	                if (i != 0)
	                    M2_L[i][j] = M2_L[i-1][j];
	                if (j != 0)
	                    M2_R[i][j] = M2_R[i][j-1];
	                continue;
		    }
	
	            double penalty = 0.0;
	            double x = peaks[i]+peaks[N-j];
		    double y = peaks2[common2[i]]+peaks2[common2[N-j]];
		    if (fabs(x-(parent_mass1+18)) < mass_tolerance ||
			fabs(y-(parent_mass2+18)) < mass_tolerance) 
			penalty = same_vertex_penalty;
	
		    if (i > j) {
			// s = 0
	                double prev_score = suffix2[j];
			for (vector<int>::const_iterator it = left_jumps[i].begin();
			     it != left_jumps[i].end(); ++it) {
			    prev_score = max(prev_score, D1[*it][j]);
			}
	                int i2 = prev[i];
	                if (i2 != None)
	                    prev_score = max(prev_score, M1_L[i2][j]);
			for (vector<int>::const_iterator it = left_jumps2[i].begin();
			     it != left_jumps2[i].end(); ++it) {
			    prev_score = max(prev_score, D2max[*it][j]);
			}
	                i2 = prev2[i];
	                if (i2 != None)
	                    prev_score = max(prev_score, M2_L[i2][j]);
	                D2[i][j][0] = common2_scores[i] + prev_score + penalty;
	
	                // s > 0
			for (int s = 0; s < l; ++s) {
			    int i2 = left_neighbors[i][s];
			    if (common2[i2] == None)
	                        continue;
	                    double penalty2 = penalty;
	                    double y = peaks2[common2[i]]+peaks2[common2[i2]];
	                    if (fabs(y-(parent_mass2+18)) < mass_tolerance)
	                        penalty2 += same_vertex_penalty;
	                    prev_score = D2max[i2][j];
	                    D2[i][j][s+1] = common2_scores[i] + prev_score + penalty2;
			}
		    } else { // i <= j
			// s = 0
	                double next_score = prefix2[i][0];
			for (vector<int>::const_iterator it = right_jumps[N-j].begin();
			     it != right_jumps[N-j].end(); ++it) {
			    next_score = max(next_score, D2[i][N-*it][0]);
			}
	                int j2 = next[N-j];
	                if (j2 != None) {
	                    next_score = max(next_score, M2_R[i][N-j2][0]);
			}
	                D2[i][j][0] = common2_scores[N-j] + next_score + penalty;
	
	                // s > 0
			for (int s = 0; s < l; ++s) {
			    int i2 = left_neighbors[i][s];
			    if (common2[i2] == None)
	                        continue;
	                    double penalty2 = penalty;
	                    double y = peaks2[common2[i2]]+peaks2[common2[N-j]];
	                    if (fabs(y-(parent_mass2+18)) < mass_tolerance)
	                        penalty2 += same_vertex_penalty;
	                    double next_score = prefix2[i][s+1];
			    for (vector<int>::const_iterator it = right_jumps[N-j].begin();
				 it != right_jumps[N-j].end(); ++it) {
				next_score = max(next_score, D2[i][N-*it][s+1]);
			    }
			    int j2 = next[N-j];
	                    if (j2 != None)
	                        next_score = max(next_score, M2_R[i][N-j2][s+1]);
	                    D2[i][j][s+1] = common2_scores[N-j] + next_score + penalty2;
			}
		    }
	
		    // Compute D2max
	            D2max[i][j] = *max_element(D2[i][j].begin(), D2[i][j].end());
	
		    // Compute M2_L
	            if (i != 0)
	                M2_L[i][j] = max(M2_L[i-1][j], D2max[i][j]);
	            else
	                M2_L[i][j] = D2max[i][j];
	
		    // compute M2_R
	            if (j != 0) {
			for (int s = 0; s < l+1; ++s) {
	                    M2_R[i][j][s] = max(M2_R[i][j-1][s], D2[i][j][s]);
			}
		    } else
	                M2_R[i][j] = D2[i][j];
		}
	    }
	
	    // compute D3/M3_R
	    for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
	            int l = right_neighbors[N-j].size();
	            D3[i][j].resize(l+1, -infinity);
	            M3_L[i][j].resize(l+1, -infinity);
		}
		if (common[i] == None) {
	            if (i > 0)
	                M3_L[i] = M3_L[i-1];
		    continue;
		}
	
	        for (int j = 0; j < n; ++j) {
	            if (common[N-j] == None) {
	                if (i != 0)
	                    M3_L[i][j] = M3_L[i-1][j];
	                if (j != 0)
	                    M3_R[i][j] = M3_R[i][j-1];
	                continue;
		    }
	
		    int l = right_neighbors[N-j].size();
	            double penalty = 0.0;
	            double x = peaks[i]+peaks[N-j];
		    double y = peaks2[common[i]]+peaks2[common[N-j]];
		    if (fabs(x-(parent_mass1+18)) < mass_tolerance ||
			fabs(y-(parent_mass2+18)) < mass_tolerance) 
			penalty = same_vertex_penalty;
	
		    if (i >= j) {
			// s = 0
	                double prev_score = suffix1[j][0];
			for (vector<int>::const_iterator it = left_jumps[i].begin();
			     it != left_jumps[i].end(); ++it) {
			    prev_score = max(prev_score, D3[*it][j][0]);
			}
	                int i2 = prev[i];
	                if (i2 != None)
	                    prev_score = max(prev_score, M3_L[i2][j][0]);
	                D3[i][j][0] = common_scores[i] + prev_score + penalty;
	                // s > 0
			for (int s = 0; s < l; ++s) {
			    int j2 = right_neighbors[N-j][s];
			    if (common[j2] == None)
	                        continue;
	                    double penalty2 = penalty;
	                    double y = peaks2[common[i]]+peaks2[common[j2]];
	                    if (fabs(y-(parent_mass2+18)) < mass_tolerance)
	                        penalty2 += same_vertex_penalty;
	                    double prev_score = suffix1[j][s+1];
			    for (vector<int>::const_iterator it = left_jumps[i].begin();
				 it != left_jumps[i].end(); ++it) {
				prev_score = max(prev_score, D3[*it][j][s+1]);
			    }
			    int i2 = prev[i];
	                    if (i2 != None)
	                        prev_score = max(prev_score, M3_L[i2][j][s+1]);
	                    D3[i][j][s+1] = common_scores[i] + prev_score + penalty2;
			}
		    } else {
			// s = 0
	                double next_score = prefix1[i];
			for (vector<int>::const_iterator it = right_jumps[N-j].begin();
			     it != right_jumps[N-j].end(); ++it) {
			    next_score = max(next_score, D1[i][N-*it]);
			}
	                int j2 = next[N-j];
	                if (j2 != None) {
	                    next_score = max(next_score, M1_R[i][N-j2]);
			}
			for (vector<int>::const_iterator it = right_jumps2[N-j].begin();
			     it != right_jumps2[N-j].end(); ++it) {
			    next_score = max(next_score, D3max[i][N-*it]);
			}
	                j2 = next2[N-j];
	                if (j2 != None) {
	                    next_score = max(next_score, M3_R[i][N-j2]);
			}
	                D3[i][j][0] = common_scores[N-j] + next_score + penalty;
	                // s > 0
			for (int s = 0; s < l; ++s) {
			    int j2 = right_neighbors[N-j][s];
			    if (common[j2] == None)
	                        continue;
	                    double penalty2 = penalty;
	                    double y = peaks2[common[N-j]]+peaks2[common[j2]];
	                    if (fabs(y-(parent_mass2+18)) < mass_tolerance)
	                        penalty2 += same_vertex_penalty;
	                    double next_score = D3max[i][N-j2];
	                    D3[i][j][s+1] = common_scores[N-j] + next_score + penalty2;
			}
		    }
	
	            // Compute D3max
	            D3max[i][j] = *std::max_element(D3[i][j].begin(), D3[i][j].end());
	
		    // Compute M3_L
	            if (i != 0) {
			for (int s = 0; s < l+1; ++s) {
	                    M3_L[i][j][s] = max(M3_L[i-1][j][s], D3[i][j][s]);
			}
		    } else
	                M3_L[i][j] = D3[i][j];
	
		    // Compute M3_R
	            if (j != 0)
	                M3_R[i][j] = max(M3_R[i][j-1], D3max[i][j]);
	            else
	                M3_R[i][j] = D3max[i][j];
		}
	    }
	
	    // Find best score
	    double best_score = 0.0;
	    int best_i = None;
	    int best_j = None;
	    int best_s = None;
	    celltype best_t = INVALID;
	    for (int i = 0; i < n; ++i) {
		for (vector<int>::const_iterator it = right_jumps[i].begin();
		     it != right_jumps[i].end(); ++it) {
		    int j = N-*it;
		    if (j <= n-1) {
			if (best_score < D1[i][j]) {
			    best_score = D1[i][j];
			    best_i = i;
			    best_j = j;
			    best_s = 0;
			    best_t = cell_D1;
			}
			if (best_score < D2max[i][j]) {
			    best_score = D2max[i][j];
			    best_i = i;
			    best_j = j;
			    best_s = max_ind(D2[i][j]);
			    best_t = cell_D2;
			}
			if (best_score < D3max[i][j]) {
			    best_score = D3max[i][j];
			    best_i = i;
			    best_j = j;
			    best_s = max_ind(D3[i][j]);
			    best_t = cell_D3;
			}
		    }
		}
	
	        int j0 = next[i];
	        if (j0 != None) {
		    int j = min(n-1, N-j0);
		    if (best_score < M1_R[i][j]) {
			best_score = M1_R[i][j];
			best_i = i;
			best_j = j;
			best_s = 0;
			best_t = cell_M1_R;
		    }
	
		    double tmp = *max_element(M2_R[i][j].begin(), M2_R[i][j].end());
		    if (best_score < tmp) {
	                best_score = tmp;
	                best_i = i;
			best_j = j;
			best_s = max_ind(M2_R[i][j]);
	                best_t = cell_M2_R;
		    }
	
		    if (best_score < M3_R[i][j]) {
			best_score = M3_R[i][j];
			best_i = i;
			best_j = j;
			best_s = 0;
			best_t = cell_M3_R;
		    }
		}
	
		if (best_score < prefix1[i]) {
		    best_score = prefix1[i];
		    best_i = i;
		    best_j = 0;
		    best_s = 0;
		    best_t = cell_prefix1;
		}
	
		if (best_score < suffix1max[i]) {
		    best_score = suffix1max[i];
		    best_j = i;
		    best_i = 0;
		    best_s = max_ind(suffix1[i]);
		    best_t = cell_suffix1;
		}
	
		if (best_score < prefix2max[i]) {
		    best_score = prefix2max[i];
		    best_i = i;
		    best_j = 0;
		    best_s = max_ind(prefix2[i]);
		    best_t = cell_prefix2;
		}
	
		if (best_score < suffix2[i]) {
		    best_score = suffix2[i];
		    best_j = i;
		    best_i = 0;
		    best_s = 0;
		    best_t = cell_suffix2;
		}
	    }
	
	    return find_paths(parent_mass1, parent_mass2,
			      common_scores, common2_scores,
			      ptm_penalty,
			      prev, next,
			      prev2, next2,
			      left_jumps, right_jumps,
			      left_jumps2, right_jumps2,
			      left_neighbors, right_neighbors,
			      best_i, best_j, best_s, best_t, best_score,
			      D1, D2, D3, D2max, D3max,
			      M1_L, M1_R, M2_L, M2_R, M3_L, M3_R,
			      prefix1, suffix1, prefix2, suffix2,
			      prefix1_L, suffix1max, suffix1_R,
			      prefix2max, prefix2_L, suffix2_R);
	}
	
	//////////////////////////////////////////////////////////////////////////////
	
	pair<double, pair2> find_paths(double parent_mass1, double parent_mass2,
				       vector<double> & common_scores,
				       vector<double> & common2_scores,
				       double ptm_penalty,
				       vector<int> & prev, vector<int> & next,
				       vector<int> & prev2, vector<int> & next2,
				       vector<vector<int> > & left_jumps,
				       vector<vector<int> > & right_jumps,
				       vector<vector<int> > & left_jumps2,
				       vector<vector<int> > & right_jumps2,
				       vector<vector<int> > & left_neighbors,
				       vector<vector<int> > & right_neighbors,
				       int best_i, int best_j, int best_s,
				       celltype best_t, double best_score,
				       vector2 & D1, vector3 & D2, vector3 & D3,
				       vector2 & D2max, vector2 & D3max,
				       vector2 & M1_L, vector2 & M1_R,
				       vector2 & M2_L, vector3 & M2_R,
				       vector3 & M3_L, vector2 & M3_R,
				       vector<double> & prefix1, vector2 & suffix1,
				       vector2 & prefix2, vector<double> & suffix2,
				       vector<double> & prefix1_L,
				       vector<double> & suffix1max, vector<double> & suffix1_R,
				       vector<double> & prefix2max, vector<double> & prefix2_L,
				       vector<double> & suffix2_R)
	{
	    int n0 = prev.size();
	    int N = n0-1;
	
	    vector<int> path1;
	    vector<int> path2;
	
	    if (best_i == None)
		return pair<double, pair2>(best_score/parent_mass1, pair2(path1, path2));
	
	    int i = best_i;
	    int j = best_j;
	    int s = best_s;
	    celltype t = best_t;
	    while (1) {
		// printf("i/j/s/t %d %d %d %d\n",i,j,s,t);
		if (t == cell_prefix1) {
		    path1.insert(path1.begin(), i);
		    double prev_score = 0.0;
	            int index = 0;
	            int next_i = 0;
		    for (vector<int>::const_iterator it = left_jumps[i].begin();
			 it != left_jumps[i].end(); ++it) {
			int i2 = *it;
	                if (prefix1[i2] > prev_score) {
	                    prev_score = prefix1[i2];
	                    index = 1;
	                    next_i = i2;
			}
		    }
	            int i2 = prev[i];
	            if (i2 != None && prefix1_L[i2] > prev_score) {
	                prev_score = prefix1_L[i2];
	                index = 2;
	                next_i = i2;
		    }
	
	            i = next_i;
	            if (index == 0) {
	                break;
		    } else if (index == 1) {
			// t is unchanged
	            } else {
	                t = cell_prefix1_L;
		    }
	
		} else if (t == cell_suffix1) {
		    path1.push_back(N-j);
	            if (s == 0) {
	                double next_score = 0.0;
	                int index = 0;
	                int next_j = 0;
			for (vector<int>::const_iterator it = right_jumps[N-j].begin();
			     it != right_jumps[N-j].end(); ++it) {
			    int j2 = *it;
			    if (suffix2[N-j2]+ptm_penalty > next_score) {
				next_score = suffix2[N-j2]+ptm_penalty;
				index = 1;
				next_j = N-j2;
			    }
			}
			int j2 = next[N-j];
			if (j2 != None && suffix2_R[N-j2]+ptm_penalty > next_score) {
			    next_score = suffix2_R[N-j2]+ptm_penalty;
			    index = 2;
			    next_j = N-j2;
			}
			for (vector<int>::const_iterator it = right_jumps2[N-j].begin();
			     it != right_jumps2[N-j].end(); ++it) {
			    int j2 = *it;
	                    if (suffix1max[N-j2] > next_score) {
	                        next_score = suffix1max[N-j2];
				index = 3;
	                        next_j = N-j2;
			    }
			}
	                j2 = next2[N-j];
	                if (j2 != None && suffix1_R[N-j2] > next_score) {
	                    next_score = suffix1_R[N-j2];
	                    index = 4;
	                    next_j = N-j2;
			}
	
	                j = next_j;
	                if (index == 0) {
	                    break;
			} else if (index == 1) {
			    t = cell_suffix2;
			} else if (index == 2) {
			    t = cell_suffix2_R;
			} else if (index == 3) {
			    // t is unchanged
	                    s = max_ind(suffix1[j]);
			} else {
	                    t = cell_suffix1_R;
			}
		    } else {
			// t is unchanged
	                j = right_neighbors[N-j][s-1];
	                j = N-j;
	                s = max_ind(suffix1[j]);
		    }
	
		} else if (t == cell_prefix2) {
		    path2.insert(path2.begin(), i);
	            if (s == 0) {
	                double prev_score = 0.0;
	                int index = 0;
	                int next_i = 0;
			for (vector<int>::const_iterator it = left_jumps[i].begin();
			     it != left_jumps[i].end(); ++it) {
			    int i2 = *it;
			    if (prefix1[i2]+ptm_penalty > prev_score) {
				prev_score = prefix1[i2]+ptm_penalty;
				index = 1;
				next_i = i2;
			    }
			}
			int i2 = prev[i];
			if (i2 != None && prefix1_L[i2]+ptm_penalty > prev_score) {
			    prev_score = prefix1_L[i2]+ptm_penalty;
			    index = 2;
			    next_i = i2;
			}
			for (vector<int>::const_iterator it = left_jumps2[i].begin();
			     it != left_jumps2[i].end(); ++it) {
			    int i2 = *it;
	                    if (prefix2max[i2] > prev_score) {
				prev_score = prefix2max[i2];
				index = 3;
				next_i = i2;
			    }
			}
	                i2 = prev2[i];
	                if (i2 != None && prefix2_L[i2] > prev_score) {
	                    prev_score = prefix2_L[i2];
	                    index = 4;
	                    next_i = i2;
			}
	
	                i = next_i;
	                if (index == 0) {
	                    break;
			} else if (index == 1) {
			    t = cell_prefix1;
			} else if (index == 2) {
			    t = cell_prefix1_L; 
			} else if (index == 3) {
			    // t is unchanged
	                    s = max_ind(prefix2[i]);
	                } else {
	                    t = cell_prefix2_L;
			}
		    } else {
			// t is unchanged
	                i = left_neighbors[i][s-1];
	                s = max_ind(prefix2[i]);
		    }
	
		} else if (t == cell_suffix2) {
		    path2.push_back(N-j);
	            double next_score = 0.0;
	            int index = 0;
	            int next_j = 0;
		    for (vector<int>::const_iterator it = right_jumps[N-j].begin();
			 it != right_jumps[N-j].end(); ++it) {
			int j2 = *it;
	                if (suffix2[N-j2] > next_score) {
	                    next_score = suffix2[N-j2];
	                    index = 1;
	                    next_j = N-j2;
			}
		    }
	            int j2 = next[N-j];
	            if (j2 != None && suffix2_R[N-j2] > next_score) {
	                next_score = suffix2_R[N-j2];
	                index = 2;
	                next_j = N-j2;
		    }
	
	            j = next_j;
	            if (index == 0) {
	                break;
		    } else if (index == 1) {
			// t is unchanged
		    } else {
	                t = cell_suffix2_R;
		    }
	
	        } else if (t == cell_prefix1_L) {
	            if (prefix1_L[i] == prefix1[i]) {
	                t = cell_prefix1;
		    } else {
	                i -= 1;
		    }
	
		} else if (t == cell_prefix2_L) {
	            if (prefix2_L[i] == prefix2max[i]) {
	                t = cell_prefix2;
	                s = max_ind(prefix2[i]);
		    } else {
	                i -= 1;
		    }
	
		} else if (t == cell_suffix1_R) {
	            if (suffix1_R[j] == suffix1max[j]) {
	                t = cell_suffix1;
	                s = max_ind(suffix1[j]);
		    } else {
	                j -= 1;
		    }
	
		} else if (t == cell_suffix2_R) {
	            if (suffix2_R[j] == suffix2[j]) {
	                t = cell_suffix2;
		    } else {
	                j -= 1;
		    }
	
		} else if (t == cell_D1) {
	            if (i >= j) {
			path1.insert(path1.begin(), i);
	                double prev_score = suffix2[j]+ptm_penalty;
	                int index = 0;
	                int next_i = 0;
			for (vector<int>::const_iterator it = left_jumps[i].begin();
			     it != left_jumps[i].end(); ++it) {
			    int i2 = *it;
	                    if (D1[i2][j] > prev_score) {
	                        prev_score = D1[i2][j];
	                        index = 1;
	                        next_i = i2;
			    }
			}
	                int i2 = prev[i];
	                if (i2 != None && M1_L[i2][j] > prev_score) {
	                    prev_score = M1_L[i2][j];
	                    index = 2;
	                    next_i = i2;
			}
	
	                i = next_i;
			// note that if index=0, then the value of i is irrelevant.
	                if (index == 0) {
	                    t = cell_suffix2;
	                } else if (index == 1) {
			    // t is unchanged
	                } else {
	                    t = cell_M1_L;
			}
		    } else {
			path2.push_back(N-j);
	                double next_score = prefix1[i]+ptm_penalty;
	                int index = 0;
	                int next_j = 0;
			for (vector<int>::const_iterator it = right_jumps[N-j].begin();
			     it != right_jumps[N-j].end(); ++it) {
			    int j2 = *it;
	                    if (D1[i][N-j2] > next_score) {
	                        next_score = D1[i][N-j2];
	                        index = 1;
	                        next_j = N-j2;
			    }
			}
	                int j2 = next[N-j];
	                if (j2 != None && M1_R[i][N-j2] > next_score) {
	                    next_score = M1_R[i][N-j2];
	                    index = 2;
	                    next_j = N-j2;
			}
	
	                j = next_j;
	                if (index == 0) {
	                    t = cell_prefix1;
			} else if (index == 1) {
			    // t is unchanged
	                } else {
	                    t = cell_M1_R;
			}
		    }
	
	        } else if (t == cell_D2) {
	            if (i > j) {
			path2.insert(path2.begin(), i);
	                if (s == 0) {
	                    double prev_score = suffix2[j];
	                    int index = 0;
	                    int next_i = 0;
			    for (vector<int>::const_iterator it = left_jumps[i].begin();
				 it != left_jumps[i].end(); ++it) {
				int i2 = *it;
	                        if (D1[i2][j] > prev_score) {
	                            prev_score = D1[i2][j];
	                            index = 1;
	                            next_i = i2;
				}
			    }
	                    int i2 = prev[i];
	                    if (i2 != None && M1_L[i2][j] > prev_score) {
	                        prev_score = M1_L[i2][j];
	                        index = 2;
	                        next_i = i2;
			    }
			    for (vector<int>::const_iterator it = left_jumps2[i].begin();
				 it != left_jumps2[i].end(); ++it) { 
				int i2 = *it;
	                        if (D2max[i2][j] > prev_score) {
	                            prev_score = D2max[i2][j];
	                            index = 3;
	                            next_i = i2;
				}
			    }
	                    i2 = prev2[i];
	                    if (i2 != None && M2_L[i2][j] > prev_score) {
	                        prev_score = M2_L[i2][j];
	                        index = 4;
	                        next_i = i2;
			    }
	
	                    i = next_i;
	                    if (index == 0) {
	                        t = cell_suffix2;
	                    } else if (index == 1) {
	                        t = cell_D1;
	                    } else if ( index == 2) {
	                        t = cell_M1_L;
	                    } else if ( index == 3) {
				// t is unchanged
	                        s = max_ind(D2[i][j]);
			    } else { 
	                        t = cell_M2_L;
			    }
			} else {
			    // t is unchanged
	                    i = left_neighbors[i][s-1];
	                    s = max_ind(D2[i][j]);
			}
		    } else { // i <= j
			path2.push_back(N-j);
	                double next_score = prefix2[i][s];
	                int index = 0;
	                int next_j = 0;
			for (vector<int>::const_iterator it = right_jumps[N-j].begin();
			     it != right_jumps[N-j].end(); ++it) {
			    int j2 = *it;
	                    if (D2[i][N-j2][s] > next_score) {
	                        next_score = D2[i][N-j2][s];
	                        index = 1;
	                        next_j = N-j2;
			    }
			}
	                int j2 = next[N-j];
	                if (j2 != None && M2_R[i][N-j2][s] > next_score) {
	                    index = 2;
	                    next_j = N-j2;
			}
	
	                j = next_j;
	                if (index == 0) {
	                    t = cell_prefix2;
	                } else if (index == 1) {
			    // t is unchanged
	                } else {
	                    t = cell_M2_R;
			}
		    }
	
	        } else if (t == cell_D3) {
	            if (i >= j) {
			path1.insert(path1.begin(), i);
	                double prev_score = suffix1[j][s];
	                int index = 0;
	                int next_i = 0;
			for (vector<int>::const_iterator it = left_jumps[i].begin();
			     it != left_jumps[i].end(); ++it) { 
			    int i2 = *it;
	                    if (D3[i2][j][s] > prev_score) {
	                        prev_score = D3[i2][j][s];
	                        index = 1;
	                        next_i = i2;
			    }
			}
	                int i2 = prev[i];
	                if (i2 != None && M3_L[i2][j][s] > prev_score) {
	                    index = 2;
	                    next_i = i2;
			}
	
	                i = next_i;
	                if (index == 0) {
	                    t = cell_suffix1;
	                } else if (index == 1) {
	                    // t is unchanged
	                } else {
	                    t = cell_M3_L;
			}
		    } else { // i < j
			path1.push_back(N-j);
	                if (s == 0) {
	                    double next_score = prefix1[i];
	                    int index = 0;
	                    int next_j = 0;
			    for (vector<int>::const_iterator it = right_jumps[N-j].begin();
				 it != right_jumps[N-j].end(); ++it) {
				int j2 = *it;
	                        if (D1[i][N-j2] > next_score) {
	                            next_score = D1[i][N-j2];
	                            index = 1;
	                            next_j = N-j2;
				}
			    }
	                    int j2 = next[N-j];
	                    if (j2 != None && M1_R[i][N-j2] > next_score) {
	                        next_score = M1_R[i][N-j2];
	                        index = 2;
	                        next_j = N-j2;
			    }
			    for (vector<int>::const_iterator it = right_jumps2[N-j].begin();
				 it != right_jumps2[N-j].end(); ++it) {
				int j2 = *it;
	                        if (D3max[i][N-j2] > next_score) {
	                            next_score = D3max[i][N-j2];
	                            index = 3;
	                            next_j = N-j2;
				}
			    }
	                    j2 = next2[N-j];
	                    if (j2 != None && M3_R[i][N-j2] > next_score) {
	                        next_score = M3_R[i][N-j2];
	                        index = 4;
	                        next_j = N-j2;
			    }
	
	                    j = next_j;
	                    if (index == 0) {
	                        t = cell_prefix1;
	                    } else if (index == 1) {
	                        t = cell_D1;
			    } else if (index == 2) {
	                        t = cell_M1_R;
	                    } else if ( index == 3) {
				// t is unchanged
	                        s = max_ind(D3[i][j]);
			    } else {
	                        t = cell_M3_R;
			    }
			} else {
			    // t is unchanged
	                    j = right_neighbors[N-j][s-1];
	                    j = N-j;
	                    s = max_ind(D3[i][j]);
			}
		    }
	
	        } else if (t == cell_M1_L) {
	            if (M1_L[i][j] == D1[i][j]) {
	                t = cell_D1;
		    } else {
	                i -= 1;
		    }
	
	        } else if (t == cell_M2_L) {
	            if (M2_L[i][j] == D2max[i][j]) {
	                t = cell_D2;
	                s = max_ind(D2[i][j]);
	            } else {
	                i -= 1;
		    }
	
	        } else if (t == cell_M3_L) {
	            if (M3_L[i][j][s] == D3[i][j][s]) {
	                t = cell_D3;
	            } else {
	                i -= 1;
		    }
	
	        } else if (t == cell_M1_R) {
	            if (M1_R[i][j] == D1[i][j]) {
	                t = cell_D1;
	            } else {
	                j -= 1;
		    }
	
	        } else if (t == cell_M2_R) {
	            if (M2_R[i][j][s] == D2[i][j][s]) {
	                t = cell_D2;
	            } else {
	                j -= 1;
		    }
	
	        } else if (t == cell_M3_R) {
	            if (M3_R[i][j] == D3max[i][j]) {
	                t = cell_D3;
	                s = max_ind(D3[i][j]);
	            } else {
	                j -= 1;
		    }
		}
	    }
	
	    return pair<double, pair2>(best_score/parent_mass1, pair2(path1, path2));
	}

};

//
//  findMatchingJumps - Finds all peaks to the left and to the right of peaks[peakIdx]
//    whose peak masses correspond to a valid jump from peaks[peakIdx]. The indices
//    of these peaks are returned in leftMatches and rightMatches.
//
void findMatchingJumps(int peakIdx, vector<double> &peaks, 
						vector<float> &jumps, float peakTol,
						vector<int> &leftMatches, vector<int> &rightMatches);


float spec_align(Spectrum *spec1, Spectrum *spec2, float peakTol,
                  Spectrum *matched1, Spectrum *matched2, int maxAAJump,
                  float sameVertexPenalty, float ptmPenalty, bool forceSymmetry,
                  bool addZPMmatches) {

	const float MIN_AA_MASS = 57.0214637230;

	if(spec1->size()==0 or spec2->size()==0) {
		if(matched1) matched1->resize(0);
		if(matched2) matched1->resize(0);
		return 0;
	}
	
	AAJumps jumps(0);  float jumpsSupremum=0;
	if(maxAAJump>0) { jumps.getjumps(maxAAJump);  jumpsSupremum = jumps.masses[jumps.masses.size()-1]+2*peakTol+.00001; }

	if (spec1->parentMass > spec2->parentMass) {
		Spectrum *tmp=spec1; spec1=spec2; spec2=tmp;
		tmp=matched1; matched1=matched2; matched2=tmp;
	}

	vector<double> peaks;   peaks.reserve(spec1->size()+1);
	vector<double> peaks2(spec2->size());   for(int i=0;i<spec2->size();i++) peaks2[i]=(*spec2)[i][0];
	vector<int> common;    common.reserve(spec1->size()+1);
	vector<int> common2;   common2.reserve(spec1->size()+1);
	vector<double> common_scores;    common_scores.reserve(spec1->size()+1);
	vector<double> common2_scores;   common2_scores.reserve(spec1->size()+1);
	vector<int> prev;      prev.reserve(spec1->size());
	vector<int> next;      next.reserve(spec1->size());
	vector<int> prev2;     prev2.reserve(spec1->size());
	vector<int> next2;     next2.reserve(spec1->size());
	vector<vector<int> > left_neighbors;   left_neighbors.reserve(spec1->size());
	vector<vector<int> > right_neighbors;  right_neighbors.reserve(spec1->size());
	

	vector<TwoValues<float> > spec1sym;   // Symmetric version of spectrum1
	if (forceSymmetry) MakeSymmetric(&spec1->peakList, spec1->parentMass, peakTol, &spec1sym);
	else { spec1sym.resize(spec1->peakList.size()); for(unsigned int i=0; i<spec1sym.size(); i++) spec1sym[i]=spec1->peakList[i]; }

	// Fill in peaks, common, common2, common_scores, common2_scores	
	int idx1,       // Index in spec1sym
	    idx2=0,     // Index in spectrum 2
	    idxPeaks=0; // Index in peaks
	float pmDelta = spec2->parentMass - spec1->parentMass;
	for(idx1=0; idx1<(int)spec1sym.size(); idx1++) {
		peaks.resize(idxPeaks+1);   common.resize(idxPeaks+1);   common2.resize(idxPeaks+1);
		common_scores.resize(idxPeaks+1);    common2_scores.resize(idxPeaks+1);
		
		// Fill in common
		if(idx2>=spec2->size()) idx2=(int)spec2->size()-1;
		for(; idx2>=0 and (*spec2)[idx2][0]>=spec1sym[idx1][0]-peakTol; idx2--);  // Move below potential matches
		for(idx2=max(0,idx2); idx2<spec2->size() and (*spec2)[idx2][0]<spec1sym[idx1][0]-peakTol; idx2++);  // Find first potential match
		common_scores[idxPeaks]=0;   common[idxPeaks]=SPECNETS_EXTRA::None;
		for(; idx2<(int)spec2->size() and (*spec2)[idx2][0]<=spec1sym[idx1][0]+peakTol; idx2++)
			if ((*spec2)[idx2][1]>common_scores[idxPeaks])
				{ common_scores[idxPeaks]=(*spec2)[idx2][1]; common[idxPeaks]=idx2; }
		if (common_scores[idxPeaks]>0) common_scores[idxPeaks]+=spec1sym[idx1][1];
		
		//Fill in common2
		if(idx2>=(int)spec2->size()) idx2=(int)spec2->size()-1;
		for(; idx2>0 and (*spec2)[idx2][0]>=spec1sym[idx1][0]+pmDelta-peakTol; idx2--);
		for(; idx2<(int)spec2->size() and (*spec2)[idx2][0]<spec1sym[idx1][0]+pmDelta-peakTol; idx2++);
		common2_scores[idxPeaks]=0;   common2[idxPeaks]=SPECNETS_EXTRA::None;
		for(; idx2<(int)spec2->size() and (*spec2)[idx2][0]<=spec1sym[idx1][0]+pmDelta+peakTol; idx2++)
			if ((*spec2)[idx2][1]>common2_scores[idxPeaks])
				{ common2_scores[idxPeaks]=(*spec2)[idx2][1]; common2[idxPeaks]=idx2; }
		if (common2_scores[idxPeaks]>0) common2_scores[idxPeaks]+=spec1sym[idx1][1];
		if (common_scores[idxPeaks]>0 or common2_scores[idxPeaks]>0) {
			peaks[idxPeaks++]=spec1sym[idx1][0];
		}
	}
	peaks.resize(idxPeaks);

	// Make peaks symmetric and change the common/common2 structures accordingly
	vector<TwoValues<float> > peaksTmp(peaks.size());   for(int i=0; i<peaks.size(); i++) peaksTmp[i].set(peaks[i],common_scores[i]);
	vector<int> indices;
	if (forceSymmetry) MakeSymmetric(&peaksTmp, spec1->parentMass, peakTol, NULL, &indices);
	else {indices.resize(peaksTmp.size()); for(unsigned int i=0; i<indices.size(); i++) indices[i]=i; }

	peaks.resize(peaksTmp.size());
	common.resize(peaksTmp.size());     common_scores.resize(peaksTmp.size());
	common2.resize(peaksTmp.size());    common2_scores.resize(peaksTmp.size());
	idx1 = spec1->size()-1;    // Keeps track of the nearby peaks in spec1
	for(int i=peaksTmp.size()-1; i>=0 ; i--) {  // iterate backwards to avoid overwriting entries in common*
		peaks[i] = peaksTmp[i][0];
		if(indices[i]>=0) {
			common[i]=common[indices[i]];    common_scores[i]=common_scores[indices[i]];  // Note that indices[i]<=i
			common2[i]=common2[indices[i]];  common2_scores[i]=common2_scores[indices[i]];  
		} else {
			// Look in spec1 for closest peak with highest score
			int j; float bestScore = 0;
			for(j=idx1; j<spec1->size() & peaks[i]+peakTol>(*spec1)[j][0]; j++); // Find peaks after peaks[i]
			for(; j>=0 & peaks[i]+peakTol<(*spec1)[j][0]; j--); // Find first peak within tolerance of peaks[i]
			idx1 = j;
			for(; j>=0 & abs(peaks[i]-(*spec1)[j][0])<=peakTol; j--) bestScore=max(bestScore,(*spec1)[j][1]);  // Get score of best peak within tolerance
			common[i]=SPECNETS_EXTRA::None;    common_scores[i]=bestScore;
			common2[i]=SPECNETS_EXTRA::None;   common2_scores[i]=bestScore;  
		}
	}
	
	// Fill in prev, prev2, next, next2
	prev.resize(peaks.size());   prev2.resize(peaks.size());
	next.resize(peaks.size());   next2.resize(peaks.size());
	left_neighbors.resize(peaks.size());   right_neighbors.resize(peaks.size());
	vector<vector<int> > left_jumps(peaks.size()),
	                     left_jumps2(peaks.size()),
	                     right_jumps(peaks.size()),
	                     right_jumps2(peaks.size());
	int idxPrev=0,      // Keeps track of the righmost peak that trails the current by >=57 Da
	    idxPrevDelta=0, // Keeps track of the righmost peak that trails the current by >=delta Da
	    idxNext=0,      // Keeps track of the leftmost peak ahead of current by >=57 Da
	    idxNextDelta=0; // Keeps track of the leftmost peak ahead of current by >=delta Da
	float delta = min(max(MIN_AA_MASS-peakTol,pmDelta+peakTol),2*MIN_AA_MASS-peakTol);
	float delta2 = jumpsSupremum-peakTol;
	for(idxPeaks=0; idxPeaks<(int)peaks.size(); idxPeaks++) {
		if(idxPrev<0) idxPrev=0;
		for(; idxPrev<peaks.size() and peaks[idxPrev]<peaks[idxPeaks]-max(MIN_AA_MASS-peakTol,delta2); idxPrev++);  idxPrev--;
		if(idxPrev>=0) prev[idxPeaks]=idxPrev; else prev[idxPeaks]=SPECNETS_EXTRA::None;

		for(; idxPrevDelta<peaks.size() and peaks[idxPrevDelta]<peaks[idxPeaks]-max(delta,delta2); idxPrevDelta++);  idxPrevDelta--;
		if(idxPrevDelta>=0) prev2[idxPeaks]=idxPrevDelta; else { prev2[idxPeaks]=SPECNETS_EXTRA::None; idxPrevDelta=0; }
		
		// Similar for next/next2
		for(idxNext=idxPeaks; idxNext<peaks.size() and peaks[idxNext]<=peaks[idxPeaks]+max(MIN_AA_MASS-peakTol,delta2); idxNext++);
		if (idxNext==peaks.size()) next[idxPeaks]=SPECNETS_EXTRA::None; else next[idxPeaks]=idxNext;
		
		for(idxNextDelta=idxNext; idxNextDelta<peaks.size() and peaks[idxNextDelta]<=peaks[idxPeaks]+max(delta,delta2); idxNextDelta++);
		if (idxNextDelta==peaks.size()) next2[idxPeaks]=SPECNETS_EXTRA::None; else next2[idxPeaks]=idxNextDelta;

		if(maxAAJump<=0) {
			if(pmDelta<=MIN_AA_MASS-2*peakTol) { left_neighbors[idxPeaks].resize(0); right_neighbors[idxPeaks].resize(0); }
			else {
				unsigned int i;  // First valid peak in the target interval
				if (prev[idxPeaks]!=SPECNETS_EXTRA::None) {
					if(prev2[idxPeaks]!=SPECNETS_EXTRA::None) left_neighbors[idxPeaks].resize(prev[idxPeaks]-prev2[idxPeaks]); else left_neighbors[idxPeaks].resize(prev[idxPeaks]+1);
					for(i=0; prev[idxPeaks]>=i and peaks[prev[idxPeaks]-i]>=peaks[idxPeaks]-delta; i++) left_neighbors[idxPeaks][i]=prev[idxPeaks]-i;
					left_neighbors[idxPeaks].resize(i);
				}
				
				// Similar for right_neighbors
				if (next[idxPeaks]!=SPECNETS_EXTRA::None) {
					if(next2[idxPeaks]!=SPECNETS_EXTRA::None) right_neighbors[idxPeaks].resize(next2[idxPeaks]-next[idxPeaks]); else right_neighbors[idxPeaks].resize(peaks.size()-next[idxPeaks]);
					for(i=0; next[idxPeaks]+i<peaks.size() and peaks[next[idxPeaks]+i]<=peaks[idxPeaks]+delta; i++) right_neighbors[idxPeaks][i]=next[idxPeaks]+i;
					right_neighbors[idxPeaks].resize(i);
				}
			}
			left_jumps[idxPeaks].resize(0);    left_jumps2[idxPeaks].resize(0);
			right_jumps[idxPeaks].resize(0);   right_jumps2[idxPeaks].resize(0);
		} else {
			findMatchingJumps(idxPeaks, peaks, jumps.masses, peakTol, left_jumps[idxPeaks], right_jumps[idxPeaks]);
			// left_neighbors[idxPeaks], right_neighbors[idxPeaks]
			unsigned int neighCount=0, jumpsCount=0; 
			left_neighbors[idxPeaks].resize(left_jumps[idxPeaks].size());   
			left_jumps2[idxPeaks].resize(left_jumps[idxPeaks].size());
			for(unsigned int i=0; i<left_jumps[idxPeaks].size(); i++) {
				if(peaks[idxPeaks]-peaks[left_jumps[idxPeaks][i]]<=delta) left_neighbors[idxPeaks][neighCount++]=left_jumps[idxPeaks][i];
				if(peaks[idxPeaks]-peaks[left_jumps[idxPeaks][i]]>=delta) left_jumps2[idxPeaks][jumpsCount++]=left_jumps[idxPeaks][i];
			}
			left_neighbors[idxPeaks].resize(neighCount);   
			left_jumps2[idxPeaks].resize(jumpsCount);

			neighCount=0, jumpsCount=0; 
			right_neighbors[idxPeaks].resize(right_jumps[idxPeaks].size());   
			right_jumps2[idxPeaks].resize(right_jumps[idxPeaks].size());
			for(unsigned int i=0; i<right_jumps[idxPeaks].size(); i++) {
				if(peaks[right_jumps[idxPeaks][i]]-peaks[idxPeaks]<=delta) right_neighbors[idxPeaks][neighCount++]=right_jumps[idxPeaks][i];
				if(peaks[right_jumps[idxPeaks][i]]-peaks[idxPeaks]>=delta) right_jumps2[idxPeaks][jumpsCount++]=right_jumps[idxPeaks][i];
			}
			right_neighbors[idxPeaks].resize(neighCount);   
			right_jumps2[idxPeaks].resize(jumpsCount);
		}
	}

	pair<double, SPECNETS_EXTRA::pair2> res;

	if(forceSymmetry)
		res = SPECNETS_EXTRA::align(spec1->parentMass-AAJumps::massMH, spec2->parentMass-AAJumps::massMH, peaks, peaks2,
				 common, common2, common_scores, common2_scores,
				 prev, next, prev2, next2, left_jumps, right_jumps, left_jumps2, right_jumps2,
				 left_neighbors, right_neighbors, sameVertexPenalty, ptmPenalty);
	else
		res = SPECNETS_EXTRA::align_simple(spec1->parentMass-AAJumps::massMH, spec2->parentMass-AAJumps::massMH, peaks, peaks2, 
				 common, common2, common_scores, common2_scores,
				 prev, next, left_jumps, right_jumps, ptmPenalty, 0);

	int idxMatch=0;
	float modPos;  // Mass value of where the mod was placed
	matched1->copyNP(*spec1);    matched2->copyNP(*spec2);
	matched1->peakList.resize(res.second.first.size()+res.second.second.size());
	matched2->peakList.resize(res.second.first.size()+res.second.second.size());

	for(unsigned int i=0; i<res.second.first.size(); i++) {
		idxPeaks = res.second.first[i];
		(*matched1)[idxMatch][0] = peaks[idxPeaks];
		if(common[idxPeaks]>=0) (*matched2)[idxMatch].set(peaks2[common[idxPeaks]], (*spec2)[common[idxPeaks]][1]);
		else (*matched2)[idxMatch].set((*matched1)[idxMatch][0], 0);
		(*matched1)[idxMatch][1] = common_scores[idxPeaks] - (*matched2)[idxMatch][1];
		idxMatch++;
	}

	for(unsigned int i=0; i<res.second.second.size(); i++) {
		idxPeaks = res.second.second[i];
		(*matched1)[idxMatch][0] = peaks[idxPeaks];
		if(common2[idxPeaks]>=0) (*matched2)[idxMatch].set(peaks2[common2[idxPeaks]], (*spec2)[common2[idxPeaks]][1]);
		else (*matched2)[idxMatch].set((*matched1)[idxMatch][0]+pmDelta, 0);
		(*matched1)[idxMatch][1] = common2_scores[idxPeaks] - (*matched2)[idxMatch][1];
		idxMatch++;
	}
	
	if(addZPMmatches and matched1->size()>0 and matched2->size()>0) {
		if(spec1->size()>1 and spec2->size()>1 and ((*matched1)[0]==(*spec1)[0] or (*matched1)[0]==(*spec1)[1]) and ((*matched2)[0]==(*spec2)[0] or (*matched2)[0]==(*spec2)[1]))
			{ matched1->peakList.resize(matched1->size()+1); for(unsigned int i=matched1->size()-1;i>0;i--) (*matched1)[i]=(*matched1)[i-1];
			  (*matched1)[0]=(*spec1)[0];  (*matched1)[1]=(*spec1)[1];
			  matched2->peakList.resize(matched2->size()+1); for(unsigned int i=matched2->size()-1;i>0;i--) (*matched2)[i]=(*matched2)[i-1];
			  (*matched2)[0]=(*spec2)[0];  (*matched2)[1]=(*spec2)[1];
			}
		if(spec1->size()>1 and spec2->size()>1 and ((*matched1)[matched1->size()-1]==(*spec1)[spec1->size()-1] or (*matched1)[matched1->size()-1]==(*spec1)[spec1->size()-2]) and 
		         ((*matched2)[matched2->size()-1]==(*spec2)[spec2->size()-1] | (*matched2)[matched2->size()-1]==(*spec2)[spec2->size()-2]))
			{ matched1->peakList.resize(matched1->size()+1); (*matched1)[matched1->size()-1]=(*spec1)[spec1->size()-1];  (*matched1)[matched1->size()-2]=(*spec1)[spec1->size()-2];
			  matched2->peakList.resize(matched2->size()+1); (*matched2)[matched2->size()-1]=(*spec2)[spec2->size()-1];  (*matched2)[matched2->size()-2]=(*spec2)[spec2->size()-2];
			}
	}
	
	if(res.second.first.size()==0) modPos=0;  // No peaks in idx1 => mod at the start
	else {
		if(res.second.second.size()==0) modPos=spec1->parentMass;  // No peaks in idx2 => mod at the end
		else {
			modPos=peaks[res.second.second[0]];  // Otherwise mod was placed at the first mass of spec1 in idx2
		}
	}
	return modPos;
}

void findMatchingJumps(int peakIdx, vector<double> &peaks, 
						vector<float> &jumps, float peakTol,
						vector<int> &leftMatches, vector<int> &rightMatches){
	int peaksPivot, jPivot, matchCount;
	float maxJump = jumps[jumps.size()-1]+peakTol;
	float peakDiff;
	
	if(peakIdx>0) {
		leftMatches.resize(peakIdx);   matchCount=0;
		for(peaksPivot=peakIdx-1, jPivot=0; peaksPivot>=0 and peaks[peakIdx]-peaks[peaksPivot]<=maxJump; peaksPivot--) {
			peakDiff = (float)peaks[peakIdx]-peaks[peaksPivot];
			// Find a jump with a mass lower than the peak difference
			for(jPivot=min(jPivot,(int)jumps.size()-1); jPivot>=0 and jumps[jPivot]>=peakDiff-peakTol; jPivot--);  jPivot=max(0,jPivot);
			// Find first jump with a mass not smaller than the peak difference - peakTol
			for(; jPivot<jumps.size() and jumps[jPivot]<peakDiff-peakTol; jPivot++);
			if(jPivot<jumps.size() and abs(jumps[jPivot]-peakDiff)<=peakTol+.00001) leftMatches[matchCount++]=peaksPivot;
		}
		leftMatches.resize(matchCount);	
	}
		
	rightMatches.resize(peaks.size()-peakIdx);   matchCount=0;
	for(peaksPivot=peakIdx+1, jPivot=0; peaksPivot<peaks.size() and peaks[peaksPivot]-peaks[peakIdx]<=maxJump; peaksPivot++) {
		peakDiff = (float)peaks[peaksPivot]-peaks[peakIdx];
		// Find a jump with a mass lower than the peak difference
		for(jPivot=min(jPivot,(int)jumps.size()-1); jPivot>=0 and jumps[jPivot]>=peakDiff-peakTol; jPivot--); jPivot=max(0,jPivot);
		// Find first jump with a mass not smaller than the peak difference - peakTol
		for(; jPivot<jumps.size() and jumps[jPivot]<peakDiff-peakTol; jPivot++);
		if(jPivot<jumps.size() and abs(jumps[jPivot]-peakDiff)<=peakTol+.00001) rightMatches[matchCount++]=peaksPivot;
	}
	rightMatches.resize(matchCount);	
}
