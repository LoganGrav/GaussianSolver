

import java.io.FileWriter;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;


public class Gaussian {
    
    public static void main(String[] args) {
    try {

        boolean sppFlag = false;    //initialize variables
        boolean dbleFlag = false;
        String fileName = null;
        String newFileName = null;
        int n = 0;

        for (String flag : args) {                      //check arguments foe filename and flags
            if (flag.equals("--spp")){
                sppFlag = true;}
            else if (flag.equals("--double")) {
                dbleFlag = true; }
            else if (flag.contains(".lin")){
                fileName = flag; 
                newFileName = fileName.replace(".lin", ".sol");
            }
        }

      
        BufferedWriter writer = new BufferedWriter(new FileWriter(newFileName));    //open file reader and writer
        BufferedReader reader = new BufferedReader(new FileReader(fileName));
        n = Integer.parseInt(reader.readLine());    //find matrix size


            double[][] dbleMatrix = new double[n][n];    //initialize matrix and vector
            double[] dbleVector = new double[n];

            float[][] fltMatrix = new float[n][n];
            float[] fltVector = new float[n];
        
        
        if (dbleFlag == true) {  //if-else used outside for loop. More code, but faster code
        for (int i = 0; i < n; i++) {
            String[] line = reader.readLine().trim().split("  ");
                for (int j = 0; j < n; j++){
                    dbleMatrix[i][j] = Double.parseDouble(line[j]);
                }
            }

            String[] line = reader.readLine().trim().split("  ");
            for (int i = 0; i < n; i++) {
                dbleVector[i] = Double.parseDouble(line[i]);
            }

          

        }

        else {
        for (int i = 0; i < n; i++) {
            String[] line = reader.readLine().trim().split("  ");
                for (int j = 0; j < n; j++){
                    fltMatrix[i][j] = Float.parseFloat(line[j]);
                }
            }

            String[] line = reader.readLine().trim().split("  ");
            for (int i = 0; i < n; i++) {
                fltVector[i] = Float.parseFloat(line[i]);
            }
            reader.close();
            writer.close();
        }
        if (sppFlag == false) {
            if (dbleFlag == true) {
                dbleGaussian(dbleMatrix, dbleVector, newFileName);
            }
            else {
                fltGaussian(fltMatrix, fltVector, newFileName);
            }
        }
        else {
            if (dbleFlag == true) {
                dbleSpp(dbleMatrix, dbleVector, newFileName);
            }
            else {
                fltSpp(fltMatrix, fltVector, newFileName);
            }
        }




    }
        catch(Exception e) {
            System.out.println("error");
        }
}

    public static void dbleGaussian(double[][] coeff, double[] con, String newFileName) {
        int n = con.length;
        double[] sol = new double[n] ;
        dbleFwdElim(coeff, con);
        dbleBackSub(coeff, con, sol);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(newFileName))) {
            for (int i = 0; i < n; i++) {
                writer.write(Double.toString(sol[i]) + " ");
        }
}
        catch (IOException e) {
            System.out.println("Output error.");
        }
        
    }

    public static void dbleFwdElim(double[][] coeff, double[] con) {
        int n = con.length;
        for (int k = 0; k < (n-1); k++) {
            for (int i = (k+1);i < n; i++) {
                double mult = (coeff[i][k])/(coeff[k][k]);
                for (int j=k; j<n; j++) {
                    coeff[i][j] = coeff[i][j] - (mult * coeff[k][j]);
                }
                con[i] = con[i] - (mult * con[k]);
            }
        }
    }


    public static void dbleBackSub(double[][] coeff, double[] con, double[] sol) {
        int n = con.length;
        sol[n-1] = con[n-1]/coeff[n-1][n-1];
        for (int i = (n-2); i >= 0; i--) {
            double sum = con[i];
            for (int j = (i+1); j < n; j++) {
                sum = sum - (coeff[i][j] * sol[j]);
            }
            sol[i] = sum/ coeff[i][i];
        }
    }







    public static void fltGaussian(float[][] coeff, float[] con, String newFileName) {
        int n = con.length;
        float[] sol = new float[n] ;
        fltFwdElim(coeff, con);
        fltBackSub(coeff, con, sol);

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(newFileName))) {
                for (int i = 0; i < n; i++) {
                    writer.write(Float.toString(sol[i]) + " ");
            }
    }
            catch (IOException e) {
                System.out.println("Output error.");
            }
    }

    public static void fltFwdElim(float[][] coeff, float[] con) {
        int n = con.length;
        for (int k = 0; k < (n-1); k++) {
            for (int i = (k+1);i < n; i++) {
                float mult = (coeff[i][k])/(coeff[k][k]);
                for (int j=k; j<n; j++) {
                    coeff[i][j] = coeff[i][j] - (mult * coeff[k][j]);
                }
                con[i] = con[i] - (mult * con[k]);
            }
        }
    }


    public static void fltBackSub(float[][] coeff, float[] con, float[] sol) {
        int n = con.length;
        sol[n-1] = con[n-1]/coeff[n-1][n-1];
        for (int i = (n-2); i >= 0; i--) {
            float sum = con[i];
            for (int j = (i+1); j < n; j++) {
                sum = sum - (coeff[i][j] * sol[j]);
            }
            sol[i] = sum/ coeff[i][i];
        }
    }










    public static void dbleSpp(double[][] coeff, double[] con, String newFileName) {
        int n = con.length;
        double[] sol = new double[n];
        int[] ind = new int[n];
            for (int i = 0; i < n; i++) {
                ind[i] = i;
            }
        dbleSppFwdElim(coeff, con, ind);
        dbleSppBackSub(coeff, con, sol, ind);

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(newFileName))) {
            for (int i = 0; i < n; i++) {
                writer.write(Double.toString(sol[i]) + " ");
        }
}
        catch (IOException e) {
            System.out.println("Output error.");
        }
    }

    public static void dbleSppFwdElim(double[][] coeff, double[] con, int[] ind) {
        int n = con.length;
        double[] scaling = new double[n];
        for (int i = 0; i < n; i++) {
            double smax = 0;
            for (int j = 1; j < n; j++) {
                if (smax < Math.abs(coeff[i][j])) {
                    smax = Math.abs(coeff[i][j]);
                }
            }
            scaling[i] = smax;
        }
        for (int k = 0; k < (n-1); k++) {
            double rmax = 0;
            int maxInd = k;
            for (int i = k; i < n; i++) {
                double r = Math.abs(coeff[ind[i]][k] / scaling[ind[i]]);
                if (r > rmax) {
                    rmax = r;
                    maxInd = i;
                }
            }
            int temp = ind[k];
            ind[k] = ind[maxInd];
            ind[maxInd] = temp;

            for (int i = (k+1); i < n; i++) {
                double mult = coeff[ind[i]][k] / coeff[ind[k]][k];
                for (int j = (k + 1); j < n; j++) {
                    coeff[ind[i]][j] = coeff[ind[i]][j] - (mult * coeff[ind[k]][j]);
                }
                con[ind[i]] = con[ind[i]] - (mult * con[ind[k]]);
            }
        }
    }







    public static void dbleSppBackSub(double[][] coeff, double[] con, double[] sol, int[] ind) {
        int n = con.length;
        double sum = 0;
            sol[n-1] = con[ind[n-1]] / coeff[ind[n-1]][n-1];
        for (int i = (n-2); i >= 0; i--) {
            sum = con[ind[i]];
            for (int j = (i + 1); j < n; j++) {
                sum = sum - (coeff[ind[i]][j] * sol[j]);
            }
            sol[i] = sum / coeff[ind[i]][i];
        }
    }



        //float spp







    public static void fltSpp(float[][] coeff, float[] con, String newFileName) {
        int n = con.length;
        float[] sol = new float[n];
        int[] ind = new int[n];
            for (int i = 0; i < n; i++) {
                ind[i] = i;
            }
        fltSppFwdElim(coeff, con, ind);
        fltSppBackSub(coeff, con, sol, ind);
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(newFileName))) {
                for (int i = 0; i < n; i++) {
                    writer.write(Float.toString(sol[i]) + " ");
            }
    }
            catch (IOException e) {
                System.out.println("Output error.");
            }
}

    public static void fltSppFwdElim(float[][] coeff, float[] con, int[] ind) {
        int n = con.length;
        float[] scaling = new float[n];
        for (int i = 0; i < n; i++) {
            float smax = 0;
            for (int j = 1; j < n; j++) {
                if (smax < Math.abs(coeff[i][j])) {
                    smax = Math.abs(coeff[i][j]);
                }
            }
            scaling[i] = smax;
        }
        for (int k = 0; k < (n-1); k++) {
            float rmax = 0;
            int maxInd = k;
            for (int i = k; i < n; i++) {
                float r = Math.abs(coeff[ind[i]][k] / scaling[ind[i]]);
                if (r > rmax) {
                    rmax = r;
                    maxInd = i;
                }
            }
            int temp = ind[k];
            ind[k] = ind[maxInd];
            ind[maxInd] = temp;

            for (int i = (k+1); i < n; i++) {
                float mult = coeff[ind[i]][k] / coeff[ind[k]][k];
                for (int j = (k + 1); j < n; j++) {
                    coeff[ind[i]][j] = coeff[ind[i]][j] - (mult * coeff[ind[k]][j]);
                }
                con[ind[i]] = con[ind[i]] - (mult * con[ind[k]]);
            }
        }
    }







    public static void fltSppBackSub(float[][] coeff, float[] con, float[] sol, int[] ind) {
        int n = con.length;
        float sum = 0;
            sol[n-1] = con[ind[n-1]] / coeff[ind[n-1]][n-1];
        for (int i = (n-2); i >= 0; i--) {
            sum = con[ind[i]];
            for (int j = (i + 1); j < n; j++) {
                sum = sum - (coeff[ind[i]][j] * sol[j]);
            }
            sol[i] = sum / coeff[ind[i]][i];
        }
    }
    }


