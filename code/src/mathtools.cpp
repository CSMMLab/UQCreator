#include "mathtools.h"

double MathTools::pythag(const double a, const double b) {
    double absa=std::fabs(a), absb=std::fabs(b);
    return (absa > absb ? absa*std::sqrt(1.0+std::pow(absb/absa,2)) :
           (absb == 0.0 ? 0.0 : absb*std::sqrt(1.0+std::pow(absa/absb,2))));
}

void MathTools::computeEigenValTriDiagMatrix(const matrix mat){
    assert( blaze::isSquare(mat) );
    int n = mat.rows();

    vector d(n),e(n);
    matrix z(n,n);

    vector d = blaze::trace(mat);
    trace(z) = 1.0;
    for(int i=0; i<n; ++i){
        i == 0 ? e(i) = 0 : e(i) = mat[i][i-1];
    }

    int m,l,iter,i,k;
    double s,r,p,g,f,dd,c,b;
    const double eps=std::numeric_limits<double>::epsilon();
    for (i=1;i<n;i++) e(i-1)=e(i);
    e(n-1)=0.0;
    for (l=0;l<n;l++) {
        iter=0;
        do {
            for (m=l;m<n-1;m++) {
                dd=std::fabs(d(m))+std::fabs(d(m+1));
                if (std::fabs(e(m)) <= eps*dd) break;
            }
            if (m != l) {
                if (iter++ == 30) throw("[computeEigenValTriDiagMatrix]: Too many iterations");
                g=(d(l+1)-d(l))/(2.0*e(l));
                r=pythag(g,1.0);
                g=d(m)-d(l)+e(l)/(g+std::copysign(r,g));
                s=c=1.0;
                p=0.0;
                for (i=m-1;i>=l;i--) {
                    f=s*e(i);
                    b=c*e(i);
                    e(i+1)=(r=pythag(f,g));
                    if (r == 0.0) {
                        d(i+1) -= p;
                        e(m)=0.0;
                        break;
                    }
                    s=f/r;
                    c=g/r;
                    g=d(i+1)-p;
                    r=(d(i)-g)*s+2.0*c*b;
                    d(i+1)=g+(p=s*r);
                    g=c*r-b;
                    for (k=0;k<n;k++) {
                        f=z(k,i+1);
                        z(k,i+1)=s*z(k,i)+c*f;
                        z(k,i)=c*z(k,i)-s*f;
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d(l) -= p;
                e(l)=g;
                e(m)=0.0;
            }
        } while (m != l);
    }


    //return std::make_pair(d,z);
}