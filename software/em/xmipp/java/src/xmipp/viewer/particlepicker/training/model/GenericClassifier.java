package xmipp.viewer.particlepicker.training.model;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Properties;

import javax.swing.JFrame;

import xmipp.jni.Classifier;
import xmipp.jni.MetaData;
import xmipp.jni.Particle;
import xmipp.utils.XmippWindowUtil;


public class GenericClassifier extends Classifier
{

	protected String classifierProperties;
	protected Properties properties;

	public GenericClassifier(String classifierProperties)
	{
		try
		{
			this.classifierProperties = classifierProperties;
			// create and load properties
			properties = new Properties();
			FileInputStream in = new FileInputStream(classifierProperties);
			properties.load(in);
			
			in.close();
			String[] paramNames = properties.getProperty("parameters").split(",");
			Parameter param;
			String label, help, value;
			params = new ArrayList<Classifier.Parameter>();
			for(String paramName : paramNames)
			{
				label = properties.getProperty(paramName + ".label");
				help = properties.getProperty(paramName + ".help");
				value = properties.getProperty(paramName + ".value");
				param = new Parameter(paramName, label, help, value);
				params.add(param);
			}
		}
		catch (IOException e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
	
	public void autopick(SupervisedPickerMicrograph micrograph)
	{
		String autopickCommand = properties.getProperty("autopickCommand");
		String convertCommand = properties.getProperty("convertCommand");
		String runDir = properties.getProperty("runDir");
		for(Classifier.Parameter param: params)
		{
			autopickCommand = autopickCommand.replace("%(" + param.name + ")", param.value);
		}
		String micPath = micrograph.getFile();
		String name = micrograph.getName();
		autopickCommand = autopickCommand.replace("%(micrograph)", micPath);
		autopickCommand = autopickCommand.replace("%(micrographName)", name);
		String output;
		try
		{
			System.out.println(autopickCommand);
			if(runDir != null)
				output = XmippWindowUtil.executeCommand(autopickCommand, true, runDir);
			else
				output = XmippWindowUtil.executeCommand(autopickCommand, true);
			System.out.println(output);
			System.out.println(convertCommand);
			output = XmippWindowUtil.executeCommand(convertCommand, true);
			System.out.println(output);
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	@Override
	public int getTrainingParticlesMinimum()
	{
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setSize(int psize)
	{
		// TODO Auto-generated method stub
		
	}


	@Override
	public boolean needsTraining()
	{
		return false;
	}
}