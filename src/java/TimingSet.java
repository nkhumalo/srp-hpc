import java.io.*;
import javax.swing.*;
import javax.swing.event.*;

public class TimingSet{

    public int numProcs;
    public int numFrames;
    public int numData;
    public double data[][][];
    public double time[];

    public TimingSet(int number, DefaultListModel list){

	JFileChooser chooser;
	JFrame dialogFrame;
	ExtensionFilter timingFilter;
    
	BufferedReader br;
	String card;

	try{

	    chooser = new JFileChooser("./");
	    timingFilter = new ExtensionFilter(".tim");
	    chooser.setFileFilter(timingFilter);
	    dialogFrame = new JFrame();
	    dialogFrame.setSize(300,400);
	    chooser.showOpenDialog(dialogFrame);

	    br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));

	    card=br.readLine();
	    numProcs  = Integer.parseInt(card.substring(1,5).trim());
	    numData   = Integer.parseInt(card.substring(6,10).trim());
	    numFrames = 0;

	    while((card=br.readLine()) != null){
		if(card.startsWith("timings")){numFrames++;};
	    };

	    br.close();

	    System.out.println("Number of frames is "+numFrames);

	    data = new double[numFrames][numProcs][numData];
	    time = new double[numFrames];

	    br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));

	    card=br.readLine();
	    for(int i=0; i<numData; i++){
		card=br.readLine();
		if(number==0) list.addElement(card);
	    };

	    for(int frame=0; frame<numFrames; frame++){
		card=br.readLine();
		card=br.readLine();
		time[frame]=Double.valueOf(card.substring(0,12)).doubleValue();
		for(int proc=0; proc<numProcs; proc++){
		    int count=0;
		    for(int i=0; i<numData; i++){
			if(count==0) card=br.readLine();
			data[frame][proc][i]=Double.valueOf(card.substring(count*7,(count+1)*7)).doubleValue();
			count++; if(count>9) count=0;
		    };
		};
	    };
	} catch(Exception e) {e.printStackTrace();};
    }

    public void plot(int index, int set, Graph graph){

	double sum;
	boolean first = true;

	for(int i=0; i<numFrames; i++){
	    sum=0.0;
	    for(int j=0; j<numProcs; j++){ sum=sum+data[i][j][index]; };
	    graph.addData(set,time[i],sum,!first,false); first=false;
	};
	graph.fillPlot();
    };

    public void sumPlot(Graph graph){

	boolean first = true;
	double[] accu = new double[numFrames];

	for(int i=0; i<numFrames; i++) accu[i]=0.0;
	for(int j=0; j<numData-1; j++){
	    first=true;
	    for(int i=0; i<numFrames; i++){
		for(int k=0; k<numProcs; k++) accu[i]=accu[i]+data[i][k][j];
		graph.addData(j,time[i],accu[i],!first,false); first=false;
	    };
	};
	graph.fillPlot();

    };

}
