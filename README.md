# Some_Ionic_Models

Here you can find a compendium of explicit and implicit schemes for the most widely-used ionic models implemented in MATLAB. I used these for different aspects of my PhD: from 2D spiral waves generation to the electromechanical simulation of the ventricles.

Its implementation in Fortran, C and Python will be later available.

The considered models are listed as follows:

- FitzHugh-Nagumo ionic model (1955-1964)

It consists in a single current obtained through  a normalised dimensionless potential \{u\} and a refractoriness variables \{v\}.

- Aliev-Panfilov (1995)

Aliev and Panfilov developed a similar phenomenological model based on the FitzHugh-Nagumo equations to mimic the electric membrane currents in a canine cardiomyocyte.

It consists in a single current obtained through  a normalised dimensionless potential \{u\} and a refractoriness variables \{v\}.

- Ten Tusscher-Noble-Noble-Panfilov (2004)

The aim of this model is to address the drawbacks of the already existing physiologic models designed at that time for the human ventricle, which were the Priebe-Beuckelmann model and its reduced version proposed by Bernus. TenTusscher et al developed an exhaustive model including recent data coming from human heart cells.

This model considers fifteen intensities, calculated by means of the electric potential \{V\}, thirteen gating variables \{m, h, j, d, f, r, s, x_{r1}, x_{r2}, x_{s}, x_{k1}, f_{Ca}, g\} and four ion concentrations \{Na, K, Ca, Ca_sr\}.

- Bueno Orovio-Cherry-Fenton (2008)

Also known as Minimal Model, it consists in a highly simplified approach created with the purpose of being applied in large realistic simulations while allowing the match of the observed behaviour in both physiologic and pathologic conditions. This model constitutes an extended version of the model proposed by Fenton and Karma where an additional variable was added with the purpose of a better recreation of the action potential shape.

This model considers three intensities, calculated by means of a normalised dimensionless potential \{u\} and three gating variables \{v, w, s\}.