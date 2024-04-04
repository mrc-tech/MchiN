# MchiN
**fiber section analysis in MATLAB**

Single MATLAB file with a class to evaluate numerically reinforced concrete sections response.

Versions:
- `MchiNstrisce`: 2D version based on section's strip discretization.


## How to use

The file `MchiNstrisce.m` (RIFERIMENTO ALL'ULTIMA VERSIONE NELLE _RELEASE_) must be located in the same folder of the MATLAB code.

```MATLAB
codice per fare le analisi di sezione, completo, commentato in maniera tale da mostrare cosa fa cosa...
```


## ToDo
Ogni volta che pubblico una release cancello questo ToDo, in maniera da tenere traccia nel README di cosa cambia tra una versione e un'altre? (non sarebbe meglio farlo con gli issue? Ma in questo modo posso tenere traccia in maniera ASCII direttamente dal mio computer).
- [ ] aggiornare il documento _LaTeX_ includendo anche le informazioni delle foto
- [ ] inventare una bella intestazione in cui metto tutte le informazioni principali dello script (forse anche il changelog all'inizio)
- [ ] eliminare i metodi e le funzioni inutili per i calcoli che deve fare la classe
- [ ] Controllare attentamente le procedure di analisi nonlineare (per la generazione della curva _momento-curvatura_ e per l'analisi del punto di stato ultimo).
- [ ] migliorare l'implementazione dei materiali (due definizioni delle equazioni nonlineari e dei parametri separatamente non e' per niente modulare)
- [ ] implementare la procedura `integrateStrain` in cui la discretizzazione in strisce cambia al variare dello stato tensionale nella sezione. _(Questa feature gia' fa riferimento alla versione 2.0 dello script, e' una modifica sostanziale)_
- [ ] nella funzione `initgeo` va inizializzata la sezione **OMOGENEIZZATA**
- [ ] modificare `yieldingPoint`
- [ ] la curvatura in `curvaMchi` è negativa se comprime le fibre superiori _(modifica per la versione 2.0)_
