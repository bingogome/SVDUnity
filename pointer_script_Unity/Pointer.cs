using System.Collections;
using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using UnityEngine.Windows.Speech;

namespace AROF
{
    public class Pointer : MonoBehaviour
    {
        private Dictionary<string, Action> keywordActions = new Dictionary<string, Action>();
        private KeywordRecognizer keywordRecognizer;
        private Vector3[] pointCloudDig = new Vector3[7];
        private Vector3[] pointCloudModel = new Vector3[7];
        int currDigIndex;

        // Start is called before the first frame update
        void Start()
        {
            keywordActions.Add("digitize", Digitize);
            keywordActions.Add("register", Register);

            keywordRecognizer = new KeywordRecognizer(keywordActions.Keys.ToArray());
            keywordRecognizer.OnPhraseRecognized += OnKeywordsRecognized;
            keywordRecognizer.Start();

            pointCloudModel[0] = new Vector3 { };
            pointCloudModel[1] = new Vector3 { };
            pointCloudModel[2] = new Vector3 { };
            pointCloudModel[3] = new Vector3 { };
            pointCloudModel[4] = new Vector3 { };
            pointCloudModel[5] = new Vector3 { };
            pointCloudModel[6] = new Vector3 { };

            currDigIndex = 0;

            // test svd
            //float[][] A = GetZero();
            //A[0][0] = 3.3f; A[0][1] = 4.2f; A[0][2] = 1f;
            //A[1][0] = 1.5f; A[1][1] = 2.0f; A[1][2] = 6.8f;
            //A[2][0] = 1.7f; A[2][1] = 11.0f; A[2][2] = 1f;
            //Tuple<float[][], float[][]> aaaa = QRSim3By3(A);
            //Tuple<float[][], float[][], float[][]> tulp = SVDSim3By3(A);
            //Debug.Log("" + aaaa.Item1[0][0] + " " + aaaa.Item1[0][1] + " " + aaaa.Item1[0][2]);
        }

        private void OnKeywordsRecognized(PhraseRecognizedEventArgs args)
        {
            keywordActions[args.text].Invoke();
        }

        private void Digitize()
        {
            // pointer tip position offset: (-145,-2.5,0) mm
            Vector3 tipOffset = new Vector3(-145.0f,-2.5f,0.0f);
            pointCloudDig[currDigIndex % 7] = transform.localPosition + transform.localRotation * tipOffset;

        }

        private void Register()
        {
//            n = size(A, 1);        % Length of Point Sets

//     % Calculate centroids of S and M point sets
//     a_centroid = (1 / n) * sum(A);
//            b_centroid = (1 / n) * sum(B);
//            a_centroid = a_centroid';
//    b_centroid = b_centroid';

//    % Point Deviations from centroid
//    A_tilda = A - a_centroid';
//    B_tilda = B - b_centroid';


//    % --------------------Find R that minimizes SSE --------------------%


//    % H Matrix for Singular Value Decomposition
//    H = zeros(3, 3);
//            x = 1; y = 2; z = 3;


//    % Build H
//    for i = 1:n
//        H = H + [A_tilda(i, x) * B_tilda(i, x), A_tilda(i, x) * B_tilda(i, y), A_tilda(i, x) * B_tilda(i, z);
//                 A_tilda(i, y) * B_tilda(i, x), A_tilda(i, y) * B_tilda(i, y), A_tilda(i, y) * B_tilda(i, z);
//            A_tilda(i, z) * B_tilda(i, x), A_tilda(i, z) * B_tilda(i, y), A_tilda(i, z) * B_tilda(i, z)];
//            end

//            % SVD Decomposition
//            [U, S, V] = svd(H);


//    % Check U,S,V matrices
//    format long
//    Actual = H;
//            Decomposed = U * S * V';
//    H == U * S * V';

//        % Calculate Rotation Matrix, R
//    R = V * transpose(U);


//% fprintf('The Determinant \n')
//    det(R);
//            R;


//    % Use Section IV of K. Arun et. al to change mirror and make valid matrix
//    % If this does not work, this algorithm cannot be used
//    if (det(R) < 0)
//                fprintf('Second condition reached \n')
//        V(:, 3) = -V(:, 3);
//            R = V * U';


//        if (det(R) < 0)
//                fprintf('Alogrithm failed')
//        end
//    end

//    % --------------------Find translation vector --------------------%


//    % Compute p, Translation
//    p = b_centroid - R * a_centroid;
        }

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
                for (int j = 0; j < i; j++)
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
        // Update is called once per frame
        void Update()
        {

        }
    }

}