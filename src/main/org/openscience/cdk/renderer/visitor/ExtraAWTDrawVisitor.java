/* Copyright (C) 2008 Gilleain Torrance <gilleain.torrance@gmail.com>
 *               2011 Egon Willighagen <egonw@users.sf.net>
*
*  Contact: cdk-devel@list.sourceforge.net
*
*  This program is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Lesser General Public License
*  as published by the Free Software Foundation; either version 2.1
*  of the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
*/
package org.openscience.cdk.renderer.visitor;

import java.awt.BasicStroke;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.util.Map;

import org.openscience.cdk.renderer.RendererModel;
import org.openscience.cdk.renderer.elements.ArrowElement;
import org.openscience.cdk.renderer.elements.IRenderingElement;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator.Scale;
import org.openscience.cdk.renderer.generators.ReactionSceneGenerator.ArrowHeadWidth;

/**
 * Extension of the {@link AWTDrawVisitor} that add rendering of
 * arrows.
 *
 * @cdk.module renderawtextra
 */
public class ExtraAWTDrawVisitor extends AWTDrawVisitor {
	
	public ExtraAWTDrawVisitor(Graphics2D graphics) {
		super(graphics);
	}
	
    /** {@inheritDoc} */
    public void visit(IRenderingElement element) {
        if (element instanceof ArrowElement)
            visit((ArrowElement)element);
        else
        	super.visit(element);
    }

    public void visit(ArrowElement line) {
    	RendererModel rendererModel = super.getRendererModel();
    	Map<Integer,BasicStroke> strokeMap = super.getStrokeMap();
    	Graphics2D graphics = super.getGraphics();
    	double scale = rendererModel.getParameter(
    			Scale.class).getValue();
    	Stroke savedStroke = graphics.getStroke();

    	int w = (int) (line.width * scale);
    	if (strokeMap.containsKey(w)) {
    		graphics.setStroke(strokeMap.get(w));
    	} else {
    		BasicStroke stroke = new BasicStroke(w);
    		graphics.setStroke(stroke);
    		strokeMap.put(w, stroke);
    	}

    	graphics.setColor(line.color);
    	int[] a = this.transformPoint(line.startX, line.startY);
    	int[] b = this.transformPoint(line.endX, line.endY);
    	graphics.drawLine(a[0], a[1], b[0], b[1]);
    	double aW = rendererModel.getParameter(
    			ArrowHeadWidth.class
    	).getValue() / scale;
    	if(line.direction){
    		int[] c = this.transformPoint(line.startX-aW, line.startY-aW);
    		int[] d = this.transformPoint(line.startX-aW, line.startY+aW);
    		graphics.drawLine(a[0], a[1], c[0], c[1]);
    		graphics.drawLine(a[0], a[1], d[0], d[1]);
    	}else{
    		int[] c = this.transformPoint(line.endX+aW, line.endY-aW);
    		int[] d = this.transformPoint(line.endX+aW, line.endY+aW);
    		graphics.drawLine(b[0], b[1], c[0], c[1]);
    		graphics.drawLine(b[0], b[1], d[0], d[1]);
    	}        
    	graphics.setStroke(savedStroke);
    }

}
