/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.normalization;

/**
 *
 * @author Juanjo Vega
 */
public interface iNormalizeListener {

    public void normalize(double min, double max);

    public void normalizeAuto();

    public void disableNormalization();
}
