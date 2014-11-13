/**
 * NoiseMap is a scientific computation plugin for OrbisGIS developed in order to
 * evaluate the noise impact on urban mobility plans. This model is
 * based on the French standard method NMPB2008. It includes traffic-to-noise
 * sources evaluation and sound propagation processing.
 *
 * This version is developed at French IRSTV Institute and at IFSTTAR
 * (http://www.ifsttar.fr/) as part of the Eval-PDU project, funded by the
 * French Agence Nationale de la Recherche (ANR) under contract ANR-08-VILL-0005-01.
 *
 * Noisemap is distributed under GPL 3 license. Its reference contact is Judicaël
 * Picaut <judicael.picaut@ifsttar.fr>. It is maintained by Nicolas Fortin
 * as part of the "Atelier SIG" team of the IRSTV Institute <http://www.irstv.fr/>.
 *
 * Copyright (C) 2011 IFSTTAR
 * Copyright (C) 2011-2012 IRSTV (FR CNRS 2488)
 *
 * Noisemap is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * Noisemap is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Noisemap. If not, see <http://www.gnu.org/licenses/>.
 *
 * For more information, please consult: <http://www.orbisgis.org/>
 * or contact directly:
 * info_at_ orbisgis.org
 */
package org.orbisgis.noisemap.h2;

import org.h2gis.h2spatialapi.DeterministicScalarFunction;
import org.orbisgis.noisemap.core.EvalTramwaySource;

/**
 * Return the dB(A) value corresponding to real speed in km/h (not commercial speed),tramway by hour parameters.
 * Ground category is one of theses (0:Grass,1:Rigid).Anti-vibration system can be floating panels placed under railways.
 * @author Nicolas Fortin
 */
public class BTW_EvalSource extends DeterministicScalarFunction {

    public BTW_EvalSource() {
        addProperty(PROP_REMARKS, "Return the dB(A) value corresponding to real speed in km/h (not commercial speed),tramway by hour parameters.\n" +
                "Ground category is one of theses (0:Grass,1:Rigid).Anti-vibration system can be floating panels placed under railways.\n" +
                "SELECT BTW_EvalSource(loadSpeed,tramway_count_per_hour,ground_category,has_anti_vibration_system)");
    }

    @Override
    public String getJavaStaticMethod() {
        return "evalSource";
    }

    /**
     *
     * @param speed Speed Km/h
     * @param tw_per_hour Tramway per hour
     * @param groundType Ground category [0-1]
     * @param has_anti_vibration Anti-vibration system
     * @return Noise level in dB(A)
     */
    public static double evalSource(final double speed, final double tw_per_hour, int groundType , boolean has_anti_vibration) {
        return EvalTramwaySource.evaluate(speed, tw_per_hour,
                groundType == 0 ? EvalTramwaySource.GROUND_TYPE.GRASS : EvalTramwaySource.GROUND_TYPE.RIGID,
                has_anti_vibration);
    }
}
