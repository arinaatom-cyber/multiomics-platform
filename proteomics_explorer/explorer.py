"""
Proteomics Explorer - Easy access to curated proteomics datasets
"""

import pandas as pd
import requests
import os
import tempfile
import gzip
import io
from pathlib import Path


class ProteomicsExplorer:
    
    SHEET_ID = "1M6hc3vmk1bNchMvEwXsIyyO5iq3mAzP877HTXzhzg38"
    CACHE_DIR = Path.home() / ".proteomics_explorer_cache"
    
    def __init__(self, verbose=True, use_cache=True):
        self.metadata = None
        self.current_data = None
        self.verbose = verbose
        self.use_cache = use_cache
        
        if use_cache:
            self.CACHE_DIR.mkdir(exist_ok=True)
        
        self._load_metadata()
    
    def _log(self, msg):
        if self.verbose:
            print(msg)
    
    def _load_metadata(self):
        url = f"https://docs.google.com/spreadsheets/d/{self.SHEET_ID}/export?format=csv"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            self.metadata = pd.read_csv(io.StringIO(response.text))
            self._log(f"✅ Каталог загружен: {len(self.metadata)} проектов")
        except requests.RequestException as e:
            raise ConnectionError(f"❌ Не удалось загрузить каталог: {e}")
    
    def list_projects(self, search=None, tissue=None, organism=None, limit=None):
        df = self.metadata.copy()
        
        if search:
            mask = df['Title'].str.contains(search, case=False, na=False)
            df = df[mask]
        
        if tissue:
            tissue_col = 'Tissue / Cell type (detailed)'
            if tissue_col in df.columns:
                mask = df[tissue_col].str.contains(tissue, case=False, na=False)
                df = df[mask]
        
        if organism:
            org_col = 'Organism'
            if org_col in df.columns:
                mask = df[org_col].str.contains(organism, case=False, na=False)
                df = df[mask]
        
        display_cols = ['Identifier', 'Title', 'Tissue / Cell type (detailed)', 
                       'Organism', 'TOTAL SAMPLES']
        available = [c for c in display_cols if c in df.columns]
        
        result = df[available]
        
        if limit:
            result = result.head(limit)
        
        return result.reset_index(drop=True)
    
    def get_info(self, identifier):
        row = self.metadata[self.metadata['Identifier'] == identifier]
        if row.empty:
            raise ValueError(f"❌ Проект {identifier} не найден")
        return row.iloc[0].dropna().to_dict()
    
    def _download_file(self, url, filename, size_bytes=None):
        cache_path = self.CACHE_DIR / filename
        if self.use_cache and cache_path.exists():
            self._log(f"📦 Используем кэш: {filename}")
            return cache_path
        
        self._log(f"📥 Скачиваю: {filename}")
        
        response = requests.get(url, stream=True, timeout=300)
        response.raise_for_status()
        
        save_path = cache_path if self.use_cache else Path(tempfile.mktemp())
        
        with open(save_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        return save_path
    
    def _read_file(self, filepath):
        filepath = Path(filepath)
        name = filepath.name.lower()
        
        if name.endswith('.gz'):
            with gzip.open(filepath, 'rt') as f:
                content = f.read()
            name = name[:-3]
            
            if name.endswith(('.txt', '.tsv')):
                return pd.read_csv(io.StringIO(content), sep='\t', low_memory=False)
            elif name.endswith('.csv'):
                return pd.read_csv(io.StringIO(content), low_memory=False)
        
        if name.endswith(('.txt', '.tsv')):
            return pd.read_csv(filepath, sep='\t', low_memory=False)
        elif name.endswith('.csv'):
            return pd.read_csv(filepath, low_memory=False)
        elif name.endswith(('.xlsx', '.xls')):
            return pd.read_excel(filepath)
        else:
            return pd.read_csv(filepath, sep='\t', low_memory=False)
    
    def list_files(self, identifier):
        api_url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/file/byProject?accession={identifier}"
        
        response = requests.get(api_url, timeout=30)
        if response.status_code != 200:
            raise ValueError(f"Проект {identifier} не найден в PRIDE")
        
        files = response.json()
        
        records = []
        for f in files:
            records.append({
                'fileName': f.get('fileName', ''),
                'fileType': f.get('fileType', ''),
                'sizeMB': round(f.get('fileSizeBytes', 0) / 1024 / 1024, 1),
                'downloadLink': f.get('downloadLink', '')
            })
        
        return pd.DataFrame(records)
    
    def load(self, identifier, file_pattern=None):
        self._log(f"📡 Загружаю {identifier}...")
        
        api_url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/file/byProject?accession={identifier}"
        
        response = requests.get(api_url, timeout=30)
        if response.status_code != 200:
            raise ValueError(f"Проект {identifier} не найден в PRIDE")
        
        files = response.json()
        
        table_exts = ('.txt', '.csv', '.tsv', '.xlsx', '.txt.gz', '.tsv.gz', '.csv.gz')
        candidates = [
            f for f in files 
            if str(f.get('fileName', '')).lower().endswith(table_exts)
            and f.get('fileSizeBytes', 0) > 50000
        ]
        
        if not candidates:
            raise FileNotFoundError(f"В {identifier} нет готовых таблиц")
        
        if file_pattern:
            pattern_match = [
                f for f in candidates 
                if file_pattern.lower() in f['fileName'].lower()
            ]
            if pattern_match:
                candidates = pattern_match
        
        priority_keywords = ['protein', 'abundance', 'intensity', 'quant', 
                           'report', 'result', 'maxquant', 'diann']
        
        best = None
        for kw in priority_keywords:
            for f in candidates:
                if kw in f['fileName'].lower():
                    best = f
                    break
            if best:
                break
        
        if not best:
            best = max(candidates, key=lambda x: x['fileSizeBytes'])
        
        filename = best['fileName']
        url = best['downloadLink']
        size = best['fileSizeBytes']
        
        self._log(f"📄 Файл: {filename} ({size/1024/1024:.1f} MB)")
        
        filepath = self._download_file(url, filename, size)
        
        self._log("📊 Читаю данные...")
        data = self._read_file(filepath)
        
        if not self.use_cache and filepath.exists():
            filepath.unlink()
        
        self._log(f"✅ Готово: {data.shape[0]:,} строк × {data.shape[1]} колонок")
        self.current_data = data
        
        return data
    
    def tissues(self):
        col = 'Tissue / Cell type (detailed)'
        if col in self.metadata.columns:
            return sorted(self.metadata[col].dropna().unique().tolist())
        return []
    
    def organisms(self):
        col = 'Organism'
        if col in self.metadata.columns:
            return sorted(self.metadata[col].dropna().unique().tolist())
        return []
    
    def clear_cache(self):
        if self.CACHE_DIR.exists():
            import shutil
            shutil.rmtree(self.CACHE_DIR)
            self.CACHE_DIR.mkdir(exist_ok=True)
            self._log("🧹 Кэш очищен")


def launch():
    try:
        import ipywidgets as widgets
        from IPython.display import display, clear_output, HTML
    except ImportError:
        raise ImportError("Установите: pip install ipywidgets")
    
    explorer = ProteomicsExplorer()
    
    output_info = widgets.Output()
    output_table = widgets.Output()
    
    options = []
    for _, row in explorer.metadata.iterrows():
        pid = row.get('Identifier')
        if pd.notna(pid) and pid:
            title = str(row.get('Title', ''))[:50]
            samples = row.get('TOTAL SAMPLES', '?')
            options.append((f"{pid} | {title}... ({samples})", pid))
    
    if not options:
        display(HTML("<p style='color:red'>❌ Каталог пуст</p>"))
        return explorer
    
    dropdown = widgets.Dropdown(
        options=options,
        description='',
        layout={'width': '700px'}
    )
    
    def on_select(change):
        with output_info:
            clear_output()
            try:
                info = explorer.get_info(change['new'])
                html = f"""
                <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                            padding: 20px; border-radius: 15px; color: white; margin: 10px 0;">
                    <h4 style="margin:0;">📊 {info.get('Title', 'N/A')}</h4>
                    <hr style="border-color: rgba(255,255,255,0.3);">
                    <b>🧬 Ткань:</b> {info.get('Tissue / Cell type (detailed)', 'N/A')}<br>
                    <b>🐁 Организм:</b> {info.get('Organism', 'N/A')}<br>
                    <b>👥 Образцы:</b> {info.get('TOTAL SAMPLES', 'N/A')}
                </div>
                """
                display(HTML(html))
            except Exception as e:
                display(HTML(f"<p style='color:red'>❌ {e}</p>"))
    
    def on_load(b):
        with output_table:
            clear_output()
            display(HTML("<p>⏳ Загрузка...</p>"))
            try:
                clear_output()
                data = explorer.load(dropdown.value)
                display(HTML("<h4>📋 Данные:</h4>"))
                display(data.head(10))
                import builtins
                builtins.df = data
                display(HTML("<p style='color:green'>✅ Сохранено в df</p>"))
            except Exception as e:
                clear_output()
                display(HTML(f"<p style='color:red'>❌ {e}</p>"))
    
    dropdown.observe(on_select, names='value')
    
    btn = widgets.Button(description="📥 ЗАГРУЗИТЬ", button_style='success')
    btn.on_click(on_load)
    
    display(HTML("<h2>🧬 Proteomics Explorer</h2>"))
    display(dropdown)
    display(output_info)
    display(btn)
    display(output_table)
    
    if options:
        on_select({'new': options[0][1]})
    
    return explorer
