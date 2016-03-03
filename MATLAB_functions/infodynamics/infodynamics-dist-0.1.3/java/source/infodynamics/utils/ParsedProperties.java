package infodynamics.utils;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;
import java.util.Set;
import java.io.FileOutputStream;


public class ParsedProperties {
	private Properties properties;
	
	public ParsedProperties(Properties props) {
		properties = props;
	}
	
	public ParsedProperties(String propsFile) throws Exception {
		Properties props = new Properties();
		props.load(new FileInputStream(propsFile));
		properties = props;
	}

	public int getIntProperty(String name) throws Exception {
		String value = properties.getProperty(name);
		if (value == null) {
			throw new Exception("Property " + name + " not found");
		}
		return Integer.parseInt(value);
	}

	public long getLongProperty(String name) throws Exception {
		String value = properties.getProperty(name);
		if (value == null) {
			throw new Exception("Property " + name + " not found");
		}
		return Long.parseLong(value);
	}

	public double getDoubleProperty(String name) throws Exception {
		String value = properties.getProperty(name);
		if (value == null) {
			throw new Exception("Property " + name + " not found");
		}
		return Double.parseDouble(value);
	}

	public boolean getBooleanProperty(String name) throws Exception {
		String value = properties.getProperty(name);
		if (value == null) {
			throw new Exception("Property " + name + " not found");
		}
		return Boolean.parseBoolean(value);
	}
	
	/**
	 * Return an integer array for a comma separated integer list,
	 * or a specification of "startIndex:endIndex"
	 * @param name
	 * @return Integer array for the property
	 */
	public int[] getIntArrayProperty(String name) {
		int[] returnValues;
		String repeatString = properties.getProperty(name);
		if (repeatString.indexOf(':') >= 0) {
			// The repeats are of the format "startIndex:endIndex"
			String[] indices = repeatString.split(":");
			int startIndex = Integer.parseInt(indices[0]);
			int endIndex = Integer.parseInt(indices[1]);
			returnValues = new int[endIndex - startIndex + 1];
			for (int i = 0; i < returnValues.length; i++) {
				returnValues[i] = startIndex + i;
			}
		} else {
			// The repeats are in a comma separated format
			String[] repeatStrings = repeatString.split(",");
			returnValues = new int[repeatStrings.length];
			for (int i = 0; i < returnValues.length; i++) {
				returnValues[i] = Integer.parseInt(repeatStrings[i]);
			}
		}
		return returnValues;
	}
	
	public String getStringProperty(String name) throws Exception {
		String value = getProperty(name);
		if (value == null) {
			throw new Exception("Property " + name + " not found");
		}
		return value;
	}
	
	/**
	 * Return a string array for a comma separated string list
	 * 
	 * @param name
	 * @return String array for the property
	 */
	public String[] getStringArrayProperty(String name) {
		String repeatString = properties.getProperty(name);
		// The repeats are in a comma separated format
		String[] repeatStrings = repeatString.split(",");
		return repeatStrings;
	}

	public String getProperty(String name) {
		return properties.getProperty(name);
	}
	
	public Properties getProperties() {
		return properties;
	}
	
	public void setProperty(String name, String value) {
		properties.setProperty(name, value);
	}
	
	public Set<Object> keySet() {
		return properties.keySet();
	}
	                                                 
	public void writeToFile(String fileName,String comment) throws IOException{
		FileOutputStream writeFile = new FileOutputStream(fileName);
		properties.store(writeFile, comment);
	}
	
	public boolean containsProperty(String name) {
		return properties.containsKey(name);
	}
}
