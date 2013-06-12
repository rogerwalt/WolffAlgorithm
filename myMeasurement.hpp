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

		T mean (int i = 0) const {
			return sums_[i] / double(n_[i]);
		}

		T variance (int i = 0) const {
			T m = mean(i);
			T m2 = squared_sums_[i] / n_[i];
			return m2 - m*m;
		}

		T error (int i = 0) const {
			return sqrt( variance(i) / n_[i] );
		}

		int bins() const { return n_.size(); }
		int samples (int i = 0) const { return n_[i]; }
        
        void clear() {
            sums_.clear();
            squared_sums_.clear();
            x_.clear();
            n_.clear();
        }

	protected:
};

template <typename T> std::ostream& operator<< (std::ostream& out, const myMeasurement<T>& m) {
	int N = m.bins()-7; // the last seven bins have less than 100 measurements
	N = N>0 ? N : 0;
	out << "Result: " << m.mean() << " +- " << m.error() << std::endl;
	out << "Bins: " << N << std::endl;
	for (int i=0;i<N;i++) {
		out << "#" << i+1 << ": number = " << m.samples(i) << " mean = " << m.mean(i) << ", error = " << m.error(i) << std::endl;
	}
	return out;
}

#endif // __MEASUREMENTS_HPP