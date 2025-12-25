"""Node 8: Advanced Analysis - PCA, ML, and Confounders"""
import json
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.feature_selection import SelectFromModel
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import roc_auc_score, confusion_matrix, accuracy_score
from scipy import stats

# Setup paths
current_dir = Path(__file__).parent
root_dir = current_dir.parent
input_file = root_dir / "02_feature_engineering" / "data_cleaned.pkl"
std_file = root_dir / "03_standardization" / "data_standardized.pkl"
metadata_file = root_dir / "global_params.json"
output_file = current_dir / "advanced_results.pkl"

# Load metadata
print(f"Loading metadata from {metadata_file}...", flush=True)
with open(metadata_file, 'r') as f:
    metadata = json.load(f)
    if 'amino_acids_Conc' in metadata:
        AMINO_ACIDS = metadata['amino_acids_Conc']
    else:
        AMINO_ACIDS = metadata.get('amino_acids', [])

print("Loading data...", flush=True)
df = pd.read_pickle(input_file)
df_std = pd.read_pickle(std_file)

class AdvancedAnalyzer:
    def __init__(self, df, df_std, conc_cols):
        self.df = df
        self.df_std = df_std
        self.conc_cols = conc_cols
        
    def perform_pca(self):
        print("Performing PCA...", flush=True)
        pca = PCA(n_components=2)
        imputer = SimpleImputer(strategy='mean')
        X_imputed = imputer.fit_transform(self.df_std[self.conc_cols])
        
        # Handle any remaining NaNs (e.g. if a column was all NaNs)
        if np.isnan(X_imputed).any():
            X_imputed = np.nan_to_num(X_imputed)
            
        components = pca.fit_transform(X_imputed)
        return {
            'pc1': components[:, 0].tolist(),
            'pc2': components[:, 1].tolist(),
            'explained_variance': pca.explained_variance_ratio_.tolist(),
            'labels': self.df['type'].astype(str).tolist(),
            'status': self.df['status'].astype(str).tolist()
        }

    def compute_correlation_matrix(self):
        print("Computing Correlation Matrix...", flush=True)
        corr = self.df[self.conc_cols].corr(method='spearman')
        return {
            'matrix': corr.values.tolist(),
            'labels': corr.columns.tolist()
        }

    def analyze_clinical_correlations(self):
        print("Analyzing Clinical Correlations...", flush=True)
        results = {'EDSS': {}, 'Duration': {}}
        for col in self.conc_cols:
            # EDSS
            if 'EDSS' in self.df.columns:
                valid_edss = self.df[[col, 'EDSS']].dropna()
                if len(valid_edss) > 10:
                    r, p = stats.spearmanr(valid_edss[col], valid_edss['EDSS'])
                    results['EDSS'][col] = {'r': float(r), 'p': float(p)}
            
            # Duration
            if 'Duration' in self.df.columns:
                valid_dur = self.df[[col, 'Duration']].dropna()
                if len(valid_dur) > 10:
                    r, p = stats.spearmanr(valid_dur[col], valid_dur['Duration'])
                    results['Duration'][col] = {'r': float(r), 'p': float(p)}
        return results

    def calculate_roc_auc(self):
        print("Calculating ROC AUC scores...", flush=True)
        results = {}
        y = (self.df['status'] == 'case').astype(int)
        for col in self.conc_cols:
            X = self.df[col].fillna(self.df[col].mean()).fillna(0)
            try:
                auc = roc_auc_score(y, X)
                results[col] = float(auc)
            except:
                results[col] = 0.5
        return results

    def train_classifier(self):
        print("Training Random Forest Classifier...", flush=True)
        target_types = ['PPMS', 'SPMS', 'RRMS']
        mask = self.df['type'].isin(target_types)
        data = self.df[mask].copy()
        
        if len(data) < 20:
            return {'error': 'Insufficient data'}
            
        X = data[self.conc_cols].fillna(data[self.conc_cols].mean()).fillna(0)
        y = data['type']
        
        le = LabelEncoder()
        y_enc = le.fit_transform(y)
        
        try:
            X_train, X_test, y_train, y_test = train_test_split(
                X, y_enc, test_size=0.3, random_state=42, stratify=y_enc
            )
            
            pipeline = Pipeline([
                ('scaler', StandardScaler()),
                ('selector', SelectFromModel(RandomForestClassifier(n_estimators=100, random_state=42))),
                ('classifier', RandomForestClassifier(random_state=42, class_weight='balanced'))
            ])
            
            param_grid = {
                'selector__threshold': ['mean', 'median'],
                'classifier__n_estimators': [100, 200],
                'classifier__max_depth': [None, 10]
            }
            
            grid_search = GridSearchCV(pipeline, param_grid, cv=5, scoring='accuracy', n_jobs=1) # n_jobs=1 for safety
            grid_search.fit(X_train, y_train)
            
            best_model = grid_search.best_estimator_
            y_pred = best_model.predict(X_test)
            acc = accuracy_score(y_test, y_pred)
            cm = confusion_matrix(y_test, y_pred)
            
            # Feature importance
            selector = best_model.named_steps['selector']
            classifier = best_model.named_steps['classifier']
            selected_indices = selector.get_support(indices=True)
            selected_features = [self.conc_cols[i] for i in selected_indices]
            importances = classifier.feature_importances_
            feature_importance = dict(zip(selected_features, importances.tolist()))
            full_importance = {col: feature_importance.get(col, 0.0) for col in self.conc_cols}
            
            return {
                'accuracy': float(acc),
                'confusion_matrix': cm.tolist(),
                'classes': le.classes_.tolist(),
                'feature_importance': full_importance,
                'best_params': grid_search.best_params_
            }
        except Exception as e:
            print(f"Classification failed: {e}", flush=True)
            return {'error': str(e)}

    def analyze_confounders(self):
        print("Analyzing Confounders...", flush=True)
        results = {'Lek': {}, 'miejsce': {}}
        
        # Drug
        col_name = 'drug' if 'drug' in self.df.columns else 'Lek'
        if col_name in self.df.columns:
            for col in self.conc_cols:
                groups = [g[col].dropna() for n, g in self.df.groupby(col_name) if len(g[col].dropna()) > 3]
                if len(groups) > 1:
                    try:
                        s, p = stats.kruskal(*groups)
                        results['Lek'][col.replace('_conc', '')] = float(p)
                    except: pass

        # Place
        col_name = 'place' if 'place' in self.df.columns else 'miejsce'
        if col_name in self.df.columns:
            for col in self.conc_cols:
                groups = [g[col].dropna() for n, g in self.df.groupby(col_name) if len(g[col].dropna()) > 3]
                if len(groups) > 1:
                    try:
                        s, p = stats.kruskal(*groups)
                        results['miejsce'][col.replace('_conc', '')] = float(p)
                    except: pass
        return results

# Run analysis
analyzer = AdvancedAnalyzer(df, df_std, AMINO_ACIDS)
results = {
    'pca': analyzer.perform_pca(),
    'correlation_matrix': analyzer.compute_correlation_matrix(),
    'clinical_correlations': analyzer.analyze_clinical_correlations(),
    'roc_auc': analyzer.calculate_roc_auc(),
    'classification': analyzer.train_classifier(),
    'confounders': analyzer.analyze_confounders()
}

pd.to_pickle(results, output_file)
print(f"Saved: {output_file.name}", flush=True)
