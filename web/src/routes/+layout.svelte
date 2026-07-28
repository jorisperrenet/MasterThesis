<script lang="ts">
  import { base } from '$app/paths';
  import SiteShell from '$lib/site-kit/SiteShell.svelte';
  import '../app.css';

  let { children } = $props();
</script>

<SiteShell
  projectName="Master Thesis"
  currentProjectId="master-thesis"
  projectHref={`${base}/`}
  sourceHref="https://github.com/jorisperrenet/MasterThesis"
  brandLogo
  logoSrc={`${base}/personal-logo.svg`}
  notice="This thesis is shared for scholarly and educational use."
>
  {@render children()}
</SiteShell>
