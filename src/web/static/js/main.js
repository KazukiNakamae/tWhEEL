
// tWhEEL Report JavaScript

class tWhEELReport {
    constructor() {
        this.currentTab = 'dashboard';
        this.data = null;
        this.init();
    }
    
    async init() {
        await this.loadData();
        this.setupEventListeners();
        this.showTab('dashboard');
    }
    
    async loadData() {
        try {
            const response = await fetch('web_data/summary.json');
            this.data = await response.json();
            this.populateDashboard();
        } catch (error) {
            console.error('Failed to load data:', error);
            this.showError('Failed to load data');
        }
    }
    
    setupEventListeners() {
        // Tab switching
        document.querySelectorAll('.nav-tab').forEach(tab => {
            tab.addEventListener('click', (e) => {
                const tabId = e.target.dataset.tab;
                this.showTab(tabId);
            });
        });
        
        // Search functionality
        document.addEventListener('input', (e) => {
            if (e.target.classList.contains('search-box')) {
                this.handleSearch(e.target.value, e.target.dataset.target);
            }
        });
        
        // Filter functionality
        document.addEventListener('change', (e) => {
            if (e.target.classList.contains('filter-select')) {
                this.handleFilter(e.target.value, e.target.dataset.filter);
            }
        });
    }
    
    showTab(tabId) {
        // Tab switching
        document.querySelectorAll('.nav-tab').forEach(tab => {
            tab.classList.remove('active');
        });
        document.querySelectorAll('.tab-content').forEach(content => {
            content.classList.remove('active');
        });
        
        document.querySelector(`[data-tab="${tabId}"]`).classList.add('active');
        document.querySelector(`#${tabId}`).classList.add('active');
        
        this.currentTab = tabId;
        
        // Tab-specific initialization
        switch(tabId) {
            case 'genes':
                this.loadGeneData();
                break;
            case 'guides':
                this.loadGuideData();
                break;
            case 'go-analysis':
                this.loadGOData();
                break;
        }
    }
    
    populateDashboard() {
        if (!this.data) return;
        
        const stats = this.data.statistics;
        
        // Update statistics cards
        this.updateStatCard('protein-count', stats.gene_prediction?.protein_count || 0);
        this.updateStatCard('total-genes', stats.annotation?.total_genes || 0);
        this.updateStatCard('multi-genes', stats.annotation?.multi_genes || 0);
        this.updateStatCard('single-genes', stats.annotation?.single_genes || 0);
        this.updateStatCard('total-gos', stats.go_analysis?.total_gos || 0);
        this.updateStatCard('fully-targetable', stats.go_analysis?.fully_targetable_percent || 0);
    }
    
    updateStatCard(id, value) {
        const element = document.getElementById(id);
        if (element) {
            element.textContent = typeof value === 'number' ? value.toLocaleString() : value;
        }
    }
    
    showError(message) {
        const errorDiv = document.createElement('div');
        errorDiv.className = 'error';
        errorDiv.textContent = message;
        document.body.insertBefore(errorDiv, document.body.firstChild);
    }
    
    async loadGeneData() {
        // Gene data loading (to be implemented)
        console.log('Loading gene data...');
    }
    
    async loadGuideData() {
        // Guide RNA data loading (to be implemented)
        console.log('Loading guide data...');
    }
    
    async loadGOData() {
        // GO data loading (to be implemented)
        console.log('Loading GO data...');
    }
    
    handleSearch(query, target) {
        // Search functionality (to be implemented)
        console.log('Search:', query, 'Target:', target);
    }
    
    handleFilter(value, filter) {
        // Filter functionality (to be implemented)
        console.log('Filter:', filter, 'Value:', value);
    }
}

// Application initialization
document.addEventListener('DOMContentLoaded', () => {
    new tWhEELReport();
});
