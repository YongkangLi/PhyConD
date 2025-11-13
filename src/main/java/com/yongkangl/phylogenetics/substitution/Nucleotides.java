/*
 * Nucleotides.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package com.yongkangl.phylogenetics.substitution;

/**
 * implements DataType for nucleotides with ambiguous characters
 *
 * @version $Id: Nucleotides.java,v 1.10 2006/08/31 14:57:24 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Nucleotides {

    public static final String JC = "JC";
	public static final String F84 = "F84";
	public static final String HKY = "HKY";
	public static final String GTR = "GTR";
	
	/**
	 * Name of data type. For XML and human reading of data type.
	 */
	public static final String DESCRIPTION = "nucleotide";

    public static final int A_STATE = 0;
	public static final int C_STATE = 1;
	public static final int G_STATE = 2;
	public static final int UT_STATE = 3;

    public static final int R_STATE = 5; // A or G
    public static final int Y_STATE = 6; // C or T

	public static final int UNKNOWN_STATE = 16;
	public static final int GAP_STATE = 17;
	
	/** 
	 * A table to translate state numbers (0-17) into character codes
	 */
}
