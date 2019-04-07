import java.lang.Math;
import java.util.Random;

//http://mathworld.wolfram.com/DiscreteFourierTransform.html
//http://fourier.eng.hmc.edu/e101/lectures/Image_Processing/node6.html

public class fft
{
    public static void main() {
        double[][][] x = new double[][][] {
            {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
            {{0,0},{0,0},{70,0},{80,0},{90,0},{0,0},{0,0},{0,0}},
            {{0,0},{0,0},{90,0},{100,0},{110,0},{0,0},{0,0},{0,0}},
            {{0,0},{0,0},{110,0},{120,0},{130,0},{0,0},{0,0},{0,0}},
            {{0,0},{0,0},{130,0},{140,0},{150,0},{0,0},{0,0},{0,0}},
            {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
            {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
            {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}}
        };
        
        print(FFT2D(x,-1));
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    public static double[][] FFT(double[][] input, int sgn) {
        int N = input.length;
        double[][] twiddle = twiddle(N, sgn);
        double[][] output = FFT1(input, twiddle);
        
        for(int i=0; i<N; i++) {
            if(sgn == -1) {
                output[i][0] = output[i][0]/Math.sqrt(N);
                output[i][1] = output[i][1]/Math.sqrt(N);
            }
            
            if(sgn == 1) {
                output[i][0] = output[i][0]/Math.sqrt(N);
                output[i][1] = output[i][1]/Math.sqrt(N);
            }
        }        
        return output;
    }
    
    public static double[][][] FFT2D(double[][][] x, int sgn) {
        int N = x.length;
        
        double[][][] y1 = new double[N][N][2];
        for(int i=0; i<N; i++) {
            y1[i] = FFT(x[i],sgn);
        }
        double[][][] y2 = transpose(y1);
        
        double[][][] y3 = new double[8][8][2];
        for(int i=0; i<N; i++) {
            y3[i] = FFT(y2[i],sgn);
        }
        double[][][] y4 = transpose(y3);
        
        return y4;
    }
    
    public static double[][] copy(double[][] x) {
        double[][] y = new double[x.length][2];
        for(int i=0; i<x.length; i++) {
            assert (x[i].length == 2) : "wrong x dim";
            
            y[i][0] = x[i][0];
            y[i][1] = x[i][1];
        }
        return y;
    }

    public static double[][] rev(double[][] x) {
        double[][] y = new double[x.length][2];
        for(int i=0; i<x.length; i++) {
            assert (x[i].length == 2) : "wrong x dim";
            
            y[i][0] = x[x.length-i-1][0];
            y[i][1] = x[x.length-i-1][1];
        }
        return y;
    }
    
    public static double[][][] transpose(double[][][] x) {
        int N = x.length;
        double[][][] y = new double[N][N][2];
        for(int i=0; i<N; i++) {
            for(int j=0; j<N; j++) {
                assert (x[i][j].length == 2) : "wrong x dim";
                y[i][j][0] = x[j][i][0];
                y[i][j][1] = x[j][i][1];
            }
        }
        return y;
    }
    
    public static void bitReverse(double[][] x, int m) {
        int N = x.length;        
        for(int i=0; i<N; ++i) {        // bit reversal 
            int j=0;
            for (int k=0; k<m; ++k) {
                j=(j << 1) | (1 & (i >> k));
            }
            
            if (j < i) {
                double temp;
                
                temp = x[i][0];
                x[i][0] = x[j][0];
                x[j][0] = temp;
                
                temp = x[i][1];
                x[i][1] = x[j][1];
                x[j][1] = temp;
            }
        }
    }
    
    public static double[][] FFT4(double[][] x) {
        int N = x.length;
        double[][] y = copy(x);
        
        int m;
        for(m=0; (1<<m)<N; m++) {
        }
        
        for(int i=(m-1); i>=0; i--) {             // for log N stages
            int n = 1<<i;
            double w = -1*Math.PI/n;
            for(int k=0; k<N; k+=2*n) {      // for N components
                for(int j=0; j<n; j++) {     // for each section
                    double[] c = new double[]{Math.cos(j*w),Math.sin(j*w)};
                                        
                    double[] temp = y[k+j+n];                    
                    y[k+j+n] = sub(y[k+j], temp);
                    y[k+j+n] = mul(y[k+j+n], c);
                    
                    y[k+j] = add(y[k+j], temp);
                }
            }
        }
        bitReverse(y,m);
        return y;
    }
    
    public static double[][] FFT3(double[][] x) {
        int N = x.length;
        double[][] y = copy(x);
        
        int m;
        for(m=0; (1<<m)<N; m++) {
        }
        
        bitReverse(y,m);
        for(int i=0; i<m; i++) {             // for log N stages
            int n = 1<<i;
            double w = -1*Math.PI/n;
            for(int k=0; k<N; k+=2*n) {      // for N components
                for(int j=0; j<n; j++) {     // for each section
                    double[] c = new double[]{Math.cos(j*w),Math.sin(j*w)};
                    
                    double[] temp = mul(y[k+j+n], c);
                    
                    y[k+j+n] = sub(y[k+j], temp);                    
                    y[k+j] = add(y[k+j], temp);
                }
            }
        }
        return y;
    }
    
    public static double[][] FFT2(double[][] x, double[][] twiddle) {        
        if(x.length == 1) {
            return x;
        }
        
        else {
            assert (x.length%2==0) : x.length+" is not a power of 2";
                       
            double[][] x_first = split(x,0,2);
            double[][] x_second = split(x,1,2);
            
            double[][] x0 = add(x_first, x_second);
            double[][] x1 = mul(sub(x_first, x_second), split(twiddle, 0, 2));
            
            double[][] firstHalf = FFT2(x0, slice(twiddle, 0, 2));
            double[][] secondHalf = FFT2(x1, slice(twiddle, 0, 2));
            
            return mix(firstHalf, secondHalf);
        }
    }
    
    public static double[][] FFT1(double[][] x, double[][] twiddle) {        
        if(x.length == 1) {
            return x;
        }
        
        else {
            assert (x.length%2==0) : x.length+" is not a power of 2";
            
            double[][] x_even = FFT1(slice(x, 0, 2), slice(twiddle, 0, 2));
            double[][] x_odd = FFT1(slice(x, 1, 2), slice(twiddle, 0, 2));
            
            double[][] firstHalf = add(x_even, mul(split(twiddle, 0, 2), x_odd));
            double[][] secondHalf = add(x_even, mul(split(twiddle, 1, 2), x_odd));
            
            return concat(firstHalf, secondHalf);
        }
    }
        
    public static double[][] DFT(double[][] input, int sgn) {
        int N = input.length;
        double[][][] mat = mat(N, sgn);
        double[][] output = new double[N][2];

        for(int i=0; i<N; i++) {
            for(int j=0; j<N; j++) {
                output[i] = add(output[i], mul(mat[i][j], input[j]));
            }
            
            if(sgn == 1) {
                output[i][0] = output[i][0]/N;
                output[i][1] = output[i][1]/N;
            }
            
            if(sgn == -1) {
                output[i][0] = output[i][0];
                output[i][1] = output[i][1];
            }
        }
        return output;
    }
    
    public static double[][][] mat(int N, int sgn) {
        double[][] twiddle = twiddle(N, sgn);
        double[][][] DFT = new double[N][N][2];
        for(int j=0; j<N; j++) {
            for(int k=0; k<N; k++) {
                DFT[j][k] = twiddle[(j*k)%N];
            }
        }
        return DFT;
    }
    
    public static double[][] twiddle(int N, int sgn) {
        double[][] twiddle = new double[N][2];
        for(int i=0; i<N; i++) {
            twiddle[i][0] = Math.cos(sgn*2*Math.PI*i/N);
            twiddle[i][1] = Math.sin(sgn*2*Math.PI*i/N);
        }
        return twiddle;
    }
    
    //https://math.stackexchange.com/questions/465070/multiplying-two-complex-numbers-using-only-three-multiplications-of-real-numbers
    public static double[] mul(double[] a, double[] b) {
        assert(a.length == 2) : "a wrong dimension";
        assert(b.length == 2) : "b wrong dimension";
        
        double[] c = new double[2];

        c[0] = a[0]*b[0] - a[1]*b[1];
        c[1] = a[0]*b[1] + a[1]*b[0];
        
        return c;
    }
    
    public static double[] add(double[] a, double[] b) {
        assert(a.length == 2) : "a wrong dimension";
        assert(b.length == 2) : "b wrong dimension";
        
        double[] c = new double[2];

        c[0] = a[0]+b[0];
        c[1] = a[1]+b[1];
        
        return c;
    }
    
    public static double[] sub(double[] a, double[] b) {
        assert(a.length == 2) : "a wrong dimension";
        assert(b.length == 2) : "b wrong dimension";
        
        double[] c = new double[2];

        c[0] = a[0]-b[0];
        c[1] = a[1]-b[1];
        
        return c;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    public static double[][] mul(double[][] x, double[][] y) {
        assert(x.length == y.length) : "x,y not equal length";
        
        double[][] z = new double[x.length][2];
        for(int i=0; i<x.length; i++) {            
            z[i] = mul(x[i],y[i]);
        }
        return z;
    }
    
    public static double[][] add(double[][] x, double[][] y) {
        assert(x.length == y.length) : "x,y not equal length";
        
        double[][] z = new double[x.length][2];
        for(int i=0; i<x.length; i++) {            
            z[i] = add(x[i],y[i]);
        }
        return z;
    }
    
    public static double[][] sub(double[][] x, double[][] y) {
        assert(x.length == y.length) : "x,y not equal length";
        
        double[][] z = new double[x.length][2];
        for(int i=0; i<x.length; i++) {            
            z[i] = sub(x[i],y[i]);
        }
        return z;
    }
    
    public static double[][] concat(double[][] x, double[][] y) {
        double[][] z = new double[x.length+y.length][2];
        
        for(int i=0; i<x.length; i++) {
            assert (x[i].length == 2) : "wrong x dim";
            
            z[i][0] = x[i][0];
            z[i][1] = x[i][1];
        }
        
        for(int i=0; i<y.length; i++) {
            assert (y[i].length == 2) : "wrong y dim";
            
            z[i+x.length][0] = y[i][0];
            z[i+x.length][1] = y[i][1];
        }
        
        return z;
    }
    
    public static double[][] mix(double[][] x, double[][] y) {
        assert (x.length == y.length) : "non equal dims";
        double[][] z = new double[x.length+y.length][2];
        
        for(int i=0; i<x.length; i++) {
            assert (x[i].length == 2) : "wrong x dim";
            assert (y[i].length == 2) : "wrong y dim";
            
            z[2*i][0] = x[i][0];
            z[2*i][1] = x[i][1];
            z[2*i+1][0] = y[i][0];
            z[2*i+1][1] = y[i][1];
        }
        return z;
    }
    
    //first or second half
    public static double[][] split(double[][] x, int a, int b) {
        assert (a<b) : "a not less than b";
        assert (x.length % b == 0) : "wrong x dim";
        
        double[][] y = new double[x.length/b][2];
        for(int i=0; i<y.length; i++) {
            assert (x[i+a*x.length/b].length == 2) : "wrong x dim";
            
            y[i][0] = x[i+a*x.length/b][0];
            y[i][1] = x[i+a*x.length/b][1];
        }
        return y;
    }
    
    //evens or odds
    public static double[][] slice(double[][] x, int a, int b) {
        assert (a<b) : "a not less than b";
        assert (x.length % b == 0) : "wrong x dim";
        
        double[][] y = new double[x.length/b][2];
        for(int i=0; i<y.length; i++) {
            assert (x[a+b*i].length == 2) : "wrong x dim";
            
            y[i][0] = x[a+b*i][0];
            y[i][1] = x[a+b*i][1];
        }
        return y;
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    public static double[][] rand(int N) {
        Random r = new Random();
        double[][] rand = new double[N][2];
        for(int i=0; i<N; i++) {
            rand[i][0] = r.nextDouble();
            rand[i][1] = r.nextDouble();
        }
        return rand;
    }
    
    public static void print(double[] x) {
        assert(x.length == 2) : "x wrong dimension";
        if(x[0]>=0) {
            System.out.print(" ");
        }
        System.out.printf("%.4f", x[0]);
        if(x[1]>=0) {
            System.out.print("+");
        }
        System.out.printf("%.4f", x[1]);
        System.out.print("i");
    }
    
    public static void print(double[][][] x) {
        for(int i=0; i<x.length; i++) {
            for(int j=0; j<x[i].length; j++) {
                System.out.print("(");
                print(x[i][j]);
                System.out.print(") ");
            }
            System.out.println();
        }
        System.out.println();
    }
    
    public static void print(double[][] x) {
        for(int i=0; i<x.length; i++) {
            System.out.print("(");
            print(x[i]);
            System.out.print(") ");
        }
        System.out.println();
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //https://arxiv.org/pdf/1412.7580.pdf
    //https://jakevdp.github.io/blog/2013/08/28/understanding-the-fft/
    //http://fourier.eng.hmc.edu/e101/lectures/Image_Processing/node5.html
    //http://en.dsplib.org/content/fft_dec_in_time/fft_dec_in_time.html
    //http://en.dsplib.org/content/fft_dec_in_freq/fft_dec_in_freq.html
    public static void FFTiter() {
        int N = 16;
        int m = 4;
        
        int[] bit = new int[]{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
        
        for(int i=0; i<N; ++i) {        // bit reversal 
            int j=0;
            for (int k=0; k<m; ++k) {
                j=(j << 1) | (1 & (i >> k));
            }
            
            if (j < i) { 
                int temp = bit[i];
                bit[i] = bit[j];
                bit[j] = temp;
            }
        }
        
        for(int i=0; i<N; i++) {
            System.out.print(bit[i]+" ");
        }
        System.out.println();


        double[][] x = new double[][]{
            {0,0},{1,0},{2,0},{3,0},
            {4,0},{5,0},{6,0},{7,0},
            {8,0},{9,0},{10,0},{11,0},
            {12,0},{13,0},{14,0},{15,0}
        };
        
        double[][] y = new double[][]{
            {0,0},{8,0},{4,0},{12,0},
            {2,0},{10,0},{6,0},{14,0},
            {1,0},{9,0},{5,0},{13,0},
            {3,0},{11,0},{7,0},{15,0}
        };
        
        for(int i=0; i<m; i++) {         // for log N stages
            double n=Math.pow(2.0,(float)i);
            double w=Math.PI/n;
            
            int k=0;
            while(k<N-1) {             // for N components 
                for(int j=0; j<n; j++) {     // for each section 
                    double c = Math.cos(-j*w); 
                    double s = Math.sin(-j*w);
                    
                    int j1=k+j;

                    double tempr = x[j1+(int)n][0]*c-x[j1+(int)n][1]*s;
                    double tempi = x[j1+(int)n][1]*c+x[j1+(int)n][0]*s;

                    x[j1+(int)n][0] = x[j1][0]-tempr;
                    x[j1+(int)n][1] = x[j1][1]-tempi;

                    x[j1][0] = x[j1][0]+tempr;
                    x[j1][1] = x[j1][1]+tempi;
                }
                k+=2*n;
            }
        }
        print(x);
        print(FFT(y,-1));
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    public static void iter1(int N) {
        int m=N;
        int h=1;
        while(h<N) {
            m >>= 1;
            int s = 0;
            int r = h<<1;
            for(int i=0; i<m; i++) {
                for(int j=0; j<h; j++) {
                    int a = s+j;
                    int b = a+h;
                    
                    System.out.println(a+" "+b);
                }
                s+=r;
            }
            h=r;
        }
    }
    
    public static void iter2(int N) {
        for(int h=1; h<N; h<<=1) {
            //length of segment
            int r = 2*h;
            //number of segments
            int m = N/r;
            for(int i=0; i<m; i++) {
                //s is the  starting point
                int s = i*r;
                for(int j=0; j<h; j++) {
                    //start+j
                    int a = s+j;
                    //start+j+half
                    int b = a+h;
                    
                    System.out.println(a+" "+b);
                }
            }
        }
    }
    
    public static void iter3(int N) {
        int r=N;
        int m=1;
        while(m<N) {
            int s = 0;
            int h = r>>1;
            for(int i=0; i<m; i++) {
                for(int j=0; j<h; j++) {
                    int a = s+j;
                    int b = a+h;
                    
                    System.out.println(a+" "+b);
                }
                s+=r;
            }
            r=h;
            m<<=1;
        }
    }
    
    public static void iter4(int N) {
        for(int m=1; m<N; m<<=1) {
            //length of segment
            int r = N/m;
            //length of half-segment
            int h = r/2;
            for(int i=0; i<m; i++) {
                //s is the  starting point
                int s = i*r;
                for(int j=0; j<h; j++) {
                    //start+j
                    int a = s+j;
                    //start+j+half
                    int b = a+h;
                    
                    System.out.println(a+" "+b);
                }
            }
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //http://www.katjaas.nl/bitreversal/bitreversal.html
    public static int bitrev2(int n, int bits) {
        int nrev, N;
        int count;
        N = 1<<bits;
        count = bits-1;   // initialize the count variable
        nrev = n;
        for(n>>=1; n!=0; n>>=1) {
            nrev <<= 1;
            nrev |= n & 1;
            count--;
        }

        nrev <<= count;
        nrev &= N - 1;

        return nrev;
    }
    
    public static int bitrev(int n, int bits) {
        int i, nrev;          // nrev will store the bit-reversed pattern
        int N = 1<<bits;      // find N: shift left 1 by the number of bits
        nrev = n;   
        for(i=1; i<bits; i++) {
            n >>= 1;
            nrev <<= 1;
            nrev |= n & 1;   // give LSB of n to nrev
        }
        nrev &= N-1;         // clear all bits more significant than N-1
        
        return nrev;
    }
}