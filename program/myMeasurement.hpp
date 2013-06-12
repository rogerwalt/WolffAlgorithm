/**
 * A class to do binning analysis with measurement values.
 * By Mauro Iazzi (iazzi@itp.phys.ethz.ch)
 */

#ifndef __MEASUREMENTS_HPP
#define __MEASUREMENTS_HPP

#include <vector>
#include <iostream>

#include <cmath>

template <typename T>
class myMeasurement {
	private:
		std::vector<T> sums_;
		std::vector<T> squared_sums_;
		std::vector<T> x_;
		std::vector<int> n_;

	public:
        /** recursive */
		void add (const T &x, size_t i = 0) {
			if (i==n_.size()) {
				sums_.push_back(x);
				squared_sums_.push_back(x*x);
				x_.push_back(T());
				n_.push_back(0);
			} else {
				sums_[i] += x;
				squared_sums_[i] += x * x;
			}
			n_[i] += 1;
			if (n_[i]%2==1) {
				x_[i] = x;
			} else {
				T nx = (x + x_[i]) / 2.0;
				x_[i] = T();
				add(nx, i+1);
			}
		}

        /** not recursive but the same as above */
		void add_plain (const T &x) {
			T nx = x;
			for (size_t i=0;;i++) {
				if (i==n_.size()) {
					sums_.push_back(nx);
					squared_sums_.push_back(nx*nx);
					x_.push_back(T());
					n_.push_back(0);
				} else {
					sums_[i] += nx;
					squared_sums_[i] += nx * nx;
				}
				n_[i] += 1;
				if (n_[i]%2==1) {
					x_[i] = nx;
					break;
				} else {
					nx = (nx + x_[i]) / 2.0;
					x_[i] = T();
				}
			}
		}

        /** calculates the mean of a bin (is the same for every bin)*/
		T mean (int i = 0) const {
			return sums_[i] / double(n_[i]);
		}
        
        /** calculates the variance of a bin */
		T variance (int i = 0) const {
			T m = mean(i);
			T m2 = squared_sums_[i] / n_[i];
			return m2 - m*m;
		}

        /** calculates the error of a bin (should converge if measurements are uncorrelated) */
		T error (int i = 0) const {
			return sqrt( variance(i) / n_[i] );
		}

        /** returns number of bins */
		int bins() const { return n_.size(); }
        
        /** returns how many samples a bin contains */
		int samples (int i = 0) const { return n_[i]; }
        
        /** clears all measurement data */
        void clear() {
            sums_.clear();
            squared_sums_.clear();
            x_.clear();
            n_.clear();
        }
        
        /** calculates the autocorrelation time */
        double time (int i = 0) const {
            return (variance(i)*n_[0]/n_[i]/variance(0)-1.0)*0.5;
        }

	protected:
};

template <typename T> std::ostream& operator<< (std::ostream& out, const myMeasurement<T>& m) {
	int N = m.bins()-7; // the last seven bins have less than 100 measurements, so we drop them
	N = N>0 ? N : 0;
    
    out << "\t\t\t\t'mean':     " << m.mean() << "," << std::endl;
    out << "\t\t\t\t'stderr':   [";
    for (int i=0;i<N-1;i++) {
        out << m.mean(i) << ",";
    }
    out << m.mean(N-1) << "]," << std::endl;
    
    out << "\t\t\t\t'autocorr': [";
    for (int i=0;i<N-1;i++) {
        out << m.time(i) << ",";
    }
    out << m.mean(N-1) << "]," << std::endl;
    
//	out << "Result: " << m.mean() << " +- " << m.error() << std::endl;
//	out << "Bins: " << N << std::endl;
    
    
//	for (int i=0;i<N;i++) {
//		out << "#" << i+1 << ": number = " << m.samples(i) << " mean = " << m.mean(i) << ", error = " << m.error(i) << std::endl;
//	}
	return out;
}

#endif // __MEASUREMENTS_HPP