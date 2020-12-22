using System.IO;
using System;

class Program
{
    private float[][] GetEye()
    {
        float[][] I = new float[3][];
        I[0] = new float[3] { 1.0f, 0.0f, 0.0f };
        I[1] = new float[3] { 0.0f, 1.0f, 0.0f };
        I[2] = new float[3] { 0.0f, 0.0f, 1.0f };
        return I;
    }

    private float[][] GetZero()
    {
        float[][] I = new float[3][];
        I[0] = new float[3] { 0.0f, 0.0f, 0.0f };
        I[1] = new float[3] { 0.0f, 0.0f, 0.0f };
        I[2] = new float[3] { 0.0f, 0.0f, 0.0f };
        return I;
    }

    private float[][] GetDiag(float[] v)
    {
        float[][] D = GetZero();
        D[0][0] = v[0];
        D[1][1] = v[1];
        D[2][2] = v[2];
        return D;
    }

    private float[][] Mult(float[][] X, float[][] Y)
    {
        float[][] A = GetZero();

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    A[i][j] += X[i][k] * Y[k][j];

        return A;
    }

    private float[][] Transp(float[][] X)
    {
        float[][] A = GetZero();
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                A[j][i] = X[i][j];
        return A;
    }

    private Tuple<float[][], float[][]> QRSim3By3(float[][] A)
    {
        float[][] Q = GetZero();
        float[][] R = GetZero();

        float[][] X = GetZero();
        float[][] Y = GetZero();
        float[][] K = GetEye();
        

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                X[i][j] = A[i][j];

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Y[i][j] = A[i][j];

        for (int i = 0; i < 3; i++)
        {
            
            for (int j = 0; j < i ; j++)
            {
                

                K[j][i] = (X[0][i] * Y[0][j] + X[1][i] * Y[1][j] + X[2][i] * Y[2][j]) / (Y[0][j] * Y[0][j] + Y[1][j] * Y[1][j] + Y[2][j] * Y[2][j]);
                
                Y[0][i] = Y[0][i] - K[j][i] * Y[0][j];
                
                Y[1][i] = Y[1][i] - K[j][i] * Y[1][j];
                
                Y[2][i] = Y[2][i] - K[j][i] * Y[2][j];
                
                
            }
        }
        

        float[] vecY = new float[3] {0.0f, 0.0f, 0.0f };
        
        
        for (int i = 0; i < 3; i++)
        {
            float n = (float)Math.Sqrt(Y[0][i] * Y[0][i] + Y[1][i] * Y[1][i] + Y[2][i] * Y[2][i]);
            Q[0][i] = Y[0][i] / n;
            
            Q[1][i] = Y[1][i] / n;
            Q[2][i] = Y[2][i] / n;
            
            vecY[i] = n;
        }
        
        
        float[][] diagY = GetDiag(vecY);
        R = Mult(diagY, K);
        

        return Tuple.Create(Q, R);
    }

    private Tuple<float[][], float[][], float[][]> SVDSim3By3(float[][] A)
    {
        // Editted from the matlab code: Faiz Khan
        // https://github.com/Sable/mcbench-benchmarks/blob/master/12674-simple-svd/svdsim.m

        // If needed, change the "3" to a variable. This method does 3by3 just for simplicity
        float[][] U = new float[3][];
        float[][] S = new float[3][];
        float[][] V = new float[3][];
        float[][] Q = new float[3][];

        float tol = 2.2204E-16f * 1024; // eps floating-point relative accuracy 
        int loopmax = 300;
        int loopcount = 0;

        U = GetEye();
        V = GetEye();
        S = Transp(A);

        float err = float.MaxValue;

        while (err>tol && loopcount<loopmax)
        {
            Tuple<float[][], float[][]> tulp = QRSim3By3(Transp(S));
            Q = tulp.Item1; S = tulp.Item2;
            U = Mult(U, Q);
            tulp = QRSim3By3(Transp(S));
            Q = tulp.Item1; S = tulp.Item2;
            V = Mult(V, Q);

            // e = triu(S,1)
            float[][] e = GetZero();
            e[0][1] = S[0][1];
            e[0][2] = S[0][2];
            e[1][2] = S[1][2];

            float E = (float)Math.Sqrt(e[0][1] * e[0][1] + e[0][2] * e[0][2] + e[1][2] * e[1][2]);
            float F = (float)Math.Sqrt(S[0][0] * S[0][0] + S[1][1] * S[1][1] + S[2][2] * S[2][2]);

            if (F == 0)
                F = 1;
            err = E / F;
            loopcount++;
        }

        float[] SS = { S[0][0], S[1][1], S[2][2] };
        for(int i = 0; i < 3; i++)
        {
            float SSi = Math.Abs(SS[i]);
            S[i][i] = SSi;
            if (SS[i] < 0)
            {
                U[0][i] = -U[0][i];
                U[1][i] = -U[1][i];
                U[2][i] = -U[2][i];
            }

        }

        return Tuple.Create(U, S, V);
    }
    static void Main()
    {
         Program p = new Program();
        float[][] A = p.GetZero();
        A[0][0] = 3.3f; A[0][1] = 4.2f; A[0][2] = 1f;
        A[1][0] = 1.5f; A[1][1] = 2.0f; A[1][2] = 6.8f;
        A[2][0] = 1.7f; A[2][1] = 11.0f; A[2][2] = 1f;
        Tuple<float[][], float[][]> aaaa = p.QRSim3By3(A);
        Tuple<float[][], float[][], float[][]> tulp = p.SVDSim3By3(A);
        Console.WriteLine(""+aaaa.Item1[0][0]+" "+aaaa.Item1[0][1]+" "+aaaa.Item1[0][2]);
        Console.WriteLine(""+aaaa.Item1[1][0]+" "+aaaa.Item1[1][1]+" "+aaaa.Item1[1][2]);
        Console.WriteLine(""+aaaa.Item1[2][0]+" "+aaaa.Item1[2][1]+" "+aaaa.Item1[2][2]);
        Console.WriteLine("");
        Console.WriteLine(""+tulp.Item1[0][0]+" "+tulp.Item1[0][1]+" "+tulp.Item1[0][2]);
        Console.WriteLine(""+tulp.Item1[1][0]+" "+tulp.Item1[1][1]+" "+tulp.Item1[1][2]);
        Console.WriteLine(""+tulp.Item1[2][0]+" "+tulp.Item1[2][1]+" "+tulp.Item1[2][2]);
        
    }
}