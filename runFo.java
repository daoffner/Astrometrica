import java.io.*;

public class runFo{

   public static void main (String [] args) throws Exception
   {
      Thread findOrb= new Thread(new fo("./fo mpc.txt"));
      findOrb.start();
      Thread.sleep(25000);
      findOrb.stop();
      findOrb.join();
   }
   public static class fo implements Runnable{
   	private String command;
      public fo(String cmd)
   	{
   		command= cmd;
   	}
   	public void run()
   	{
         String s= null;
         try{
            Runtime rt = Runtime.getRuntime();
            System.out.println(rt);
            System.out.println(command);
            Process p= new ProcessBuilder(command).start();
            //Process p = rt.exec("ls");
            //Process pr = rt.exec("./fo mpc.txt");
            System.out.println(p);
            //Process pr = rt.exec("/bin/bash -c python ExcelReader.py");
            //Process pr = rt.exec("cmd /c "+command);
            //Process pr = rt.exec(new String[] { "sh",  "-c", command });
            //pr.waitFor();
            //BufferedReader buf = new BufferedReader(new InputStreamReader(pr.getInputStream()));
            //String line = "";
            //while ((line=buf.readLine())!=null) {
            //System.out.println(line);   	   
            //}
         }   
         catch(Exception e){
         
         }
      }
   }
}